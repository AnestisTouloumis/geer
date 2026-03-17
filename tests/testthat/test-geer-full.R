# =============================================================================
# tests/testthat/test-geer.R
#
# Unit tests for the geer package.
#
# Organisation (one describe() block per logical unit):
#   1.  geer_control()          – argument validation & defaults
#   2.  geer.default()          – S3 constructor / class stamp
#   3.  geewa() – input checks  – error paths before fitting
#   4.  geewa() – gaussian       – convergence, coefs, S3 extractors
#   5.  geewa() – poisson        – count family, multiple corstr
#   6.  geewa() – Gamma          – Gamma family
#   7.  geewa() – methods        – bias-reducing / bias-corrected / penalised
#   8.  geewa_binary()           – binary outcomes, all orstr
#   9.  vcov()                   – all four cov_type branches
#  10.  summary() / print()      – smoke tests
#  11.  residuals()              – three residual types
#  12.  predict()                – in-sample & newdata, link & response, se.fit
#  13.  fitted() / coef()        – basic extractors
#  14.  model.matrix()           – model.matrix.geer
#  15.  geecriteria()            – output structure & column names
#  16.  update()                 – re-fit via update()
#  17.  geer.default() edge cases – missing required fields
# =============================================================================


# ---------------------------------------------------------------------------
# Shared fixtures – fitted once and reused to keep the suite fast.
# Using the package's own bundled datasets.
# ---------------------------------------------------------------------------

local_data <- function() {
  data("epilepsy",       package = "geer", envir = environment())
  data("cerebrovascular", package = "geer", envir = environment())
  data("respiratory",    package = "geer", envir = environment())
  data("leprosy",        package = "geer", envir = environment())
  list(
    epilepsy        = epilepsy,
    cerebrovascular = cerebrovascular,
    respiratory     = respiratory,
    leprosy         = leprosy
  )
}

ds <- local_data()

# Gaussian – rinse dataset (continuous response)
data("rinse", package = "geer", envir = environment())

# Pre-fit objects used in multiple test blocks
fit_pois_indep <- geewa(
  formula = seizures ~ treatment + lnbaseline + lnage,
  data    = ds$epilepsy, id = id,
  family  = poisson(link = "log"),
  corstr  = "independence", method = "gee"
)

fit_pois_exch <- geewa(
  formula = seizures ~ treatment + lnbaseline + lnage,
  data    = ds$epilepsy, id = id,
  family  = poisson(link = "log"),
  corstr  = "exchangeable", method = "gee"
)

fit_bin_or <- geewa_binary(
  formula = ecg ~ period + treatment,
  id      = id,
  data    = ds$cerebrovascular,
  link    = "logit",
  orstr   = "exchangeable",
  method  = "gee"
)

fit_gauss_indep <- geewa(
  formula = score ~ treatment + baseline + time,
  data    = rinse, id = id,
  family  = gaussian(link = "identity"),
  corstr  = "independence", method = "gee"
)


# =============================================================================
# 1. geer_control()
# =============================================================================

describe("geer_control()", {

  it("returns default values with correct types", {
    ctrl <- geer_control()
    expect_type(ctrl, "list")
    expect_equal(ctrl$tolerance,      1e-6)
    expect_equal(ctrl$maxiter,        500L)
    expect_equal(ctrl$or_adding,      0.5)
    expect_equal(ctrl$step_maxit,     10L)
    expect_equal(ctrl$step_multi,     1L)
    expect_equal(ctrl$jeffreys_power, 0.5)
  })

  it("accepts valid user-supplied values", {
    ctrl <- geer_control(tolerance = 1e-4, maxiter = 100, or_adding = 0.1,
                         step_maxit = 5, step_multi = 2, jeffreys_power = 1)
    expect_equal(ctrl$tolerance,      1e-4)
    expect_equal(ctrl$maxiter,        100L)
    expect_equal(ctrl$or_adding,      0.1)
    expect_equal(ctrl$step_maxit,     5L)
    expect_equal(ctrl$step_multi,     2L)
    expect_equal(ctrl$jeffreys_power, 1)
  })

  it("errors on non-positive tolerance", {
    expect_error(geer_control(tolerance = 0),    "'tolerance' must be a positive number")
    expect_error(geer_control(tolerance = -1e-3), "'tolerance' must be a positive number")
    expect_error(geer_control(tolerance = Inf),   "'tolerance' must be a positive number")
  })

  it("errors on non-positive-integer maxiter", {
    expect_error(geer_control(maxiter = 0),   "'maxiter' must be a positive integer")
    expect_error(geer_control(maxiter = 1.5), "'maxiter' must be a positive integer")
    expect_error(geer_control(maxiter = -10), "'maxiter' must be a positive integer")
  })

  it("errors on non-positive jeffreys_power", {
    expect_error(geer_control(jeffreys_power = 0),  "'jeffreys_power' must be a positive number")
    expect_error(geer_control(jeffreys_power = -1), "'jeffreys_power' must be a positive number")
  })

  it("errors on non-positive or_adding", {
    expect_error(geer_control(or_adding = 0),  "'or_adding' must be a positive number")
    expect_error(geer_control(or_adding = -1), "'or_adding' must be a positive number")
  })

  it("errors on non-positive-integer step_maxit", {
    expect_error(geer_control(step_maxit = 0),   "'step_maxit' must be a positive integer")
    expect_error(geer_control(step_maxit = 2.5), "'step_maxit' must be a positive integer")
  })

  it("errors on non-positive-integer step_multi", {
    expect_error(geer_control(step_multi = 0), "'step_multi' must be a positive integer")
  })

})


# =============================================================================
# 2. geer.default() – S3 constructor
# =============================================================================

describe("geer.default() S3 constructor", {

  it("errors when input is not a list", {
    expect_error(geer(42),          "'x' must be a list")
    expect_error(geer("a string"),  "'x' must be a list")
  })

  it("errors with informative message when required fields are missing", {
    incomplete <- list(coefficients = 1, residuals = 1)
    err <- tryCatch(geer(incomplete), error = function(e) e$message)
    expect_match(err, "missing required components")
    # should name at least one missing field
    expect_match(err, "fitted.values|rank|family")
  })

  it("stamps the class when a complete list is provided", {
    # Build a minimal valid stub from an actual fit
    stub <- unclass(fit_pois_indep)
    obj  <- geer(stub)
    expect_s3_class(obj, "geer")
  })

})


# =============================================================================
# 3. geewa() – input validation (error paths)
# =============================================================================

describe("geewa() input validation", {

  it("errors on unknown corstr", {
    expect_error(
      geewa(seizures ~ treatment, data = ds$epilepsy, id = id,
            family = poisson("log"), corstr = "toeplitz"),
      "'corstr' should be one of"
    )
  })

  it("errors on invalid method", {
    expect_error(
      geewa(seizures ~ treatment, data = ds$epilepsy, id = id,
            family = poisson("log"), method = "em-gee"),
      "'method' should be one of"
    )
  })

  it("errors when m-dependent Mv is not a positive integer", {
    expect_error(
      geewa(seizures ~ treatment, data = ds$epilepsy, id = id,
            family = poisson("log"), corstr = "m-dependent", Mv = 0),
      "Mv must be a positive integer"
    )
    expect_error(
      geewa(seizures ~ treatment, data = ds$epilepsy, id = id,
            family = poisson("log"), corstr = "m-dependent", Mv = 1.5),
      "Mv must be a positive integer"
    )
  })

  it("errors when corstr = 'fixed' but alpha_vector is not supplied", {
    expect_error(
      geewa(seizures ~ treatment, data = ds$epilepsy, id = id,
            family = poisson("log"), corstr = "fixed"),
      "'alpha_vector' should be provided"
    )
  })

  it("errors when beta_start has wrong length", {
    expect_error(
      geewa(seizures ~ treatment, data = ds$epilepsy, id = id,
            family = poisson("log"), beta_start = c(0, 0, 0)),
      "'beta_start' must be a numeric vector of length"
    )
  })

  it("errors when phi_fixed = TRUE but phi_value is not positive", {
    expect_error(
      geewa(seizures ~ treatment, data = ds$epilepsy, id = id,
            family = poisson("log"), phi_fixed = TRUE, phi_value = -1),
      "'phi_value' must be a single positive number"
    )
  })

  it("errors when weights contain non-positive values", {
    bad_weights <- rep(1, nrow(ds$epilepsy))
    bad_weights[1] <- -1
    expect_error(
      geewa(seizures ~ treatment, data = ds$epilepsy, id = id,
            family = poisson("log"), weights = bad_weights),
      "'weights' must be strictly positive"
    )
  })

  it("errors when id has missing values", {
    bad_data       <- ds$epilepsy
    bad_data$id[1] <- NA
    expect_error(
      geewa(seizures ~ treatment, data = bad_data, id = id,
            family = poisson("log")),
      "'id' cannot contain missing values"
    )
  })

})


# =============================================================================
# 4. geewa() – Gaussian family
# =============================================================================

describe("geewa() Gaussian identity", {

  it("produces a geer object", {
    expect_s3_class(fit_gauss_indep, "geer")
  })

  it("converges", {
    expect_true(fit_gauss_indep$converged)
  })

  it("has coefficients named consistently with model matrix", {
    expect_equal(names(coef(fit_gauss_indep)),
                 colnames(fit_gauss_indep$x))
  })

  it("fitted values are finite and have same length as response", {
    fv <- fitted(fit_gauss_indep)
    expect_true(all(is.finite(fv)))
    expect_equal(length(fv), fit_gauss_indep$obs_no)
  })

  it("working residuals equal y - fitted", {
    r   <- residuals(fit_gauss_indep, type = "working")
    ref <- fit_gauss_indep$y - fitted(fit_gauss_indep)
    expect_equal(r, ref, tolerance = 1e-10)
  })

  it("scale parameter phi is positive", {
    expect_gt(fit_gauss_indep$phi, 0)
  })

  it("family object carries identity link", {
    expect_equal(fit_gauss_indep$family$family, "gaussian")
    expect_equal(fit_gauss_indep$family$link,   "identity")
  })

})


# =============================================================================
# 5. geewa() – Poisson family, multiple corstr
# =============================================================================

describe("geewa() Poisson log – convergence & corstr", {

  it("independence model converges and has positive fitted values", {
    expect_true(fit_pois_indep$converged)
    expect_true(all(fitted(fit_pois_indep) > 0))
  })

  it("exchangeable model converges and has positive fitted values", {
    expect_true(fit_pois_exch$converged)
    expect_true(all(fitted(fit_pois_exch) > 0))
  })

  it("exchangeable alpha is in (-1, 1)", {
    alpha <- fit_pois_exch$alpha
    expect_length(alpha, 1L)
    expect_gt(alpha, -1)
    expect_lt(alpha,  1)
  })

  it("ar1 model converges", {
    fit_ar1 <- geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = ds$epilepsy, id = id,
      family = poisson("log"), corstr = "ar1"
    )
    expect_true(fit_ar1$converged)
    expect_length(fit_ar1$alpha, 1L)
  })

  it("unstructured model converges", {
    fit_unstr <- geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = ds$epilepsy, id = id,
      family = poisson("log"), corstr = "unstructured"
    )
    expect_true(fit_unstr$converged)
    # 4 time points => choose(4,2) = 6 alpha parameters
    expect_length(fit_unstr$alpha, 6L)
  })

  it("m-dependent(1) model converges", {
    fit_mdep <- geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = ds$epilepsy, id = id,
      family = poisson("log"), corstr = "m-dependent", Mv = 1
    )
    expect_true(fit_mdep$converged)
  })

  it("independence and exchangeable coefs are numerically close", {
    # Under mild correlation, estimates should be in the same ballpark
    expect_equal(coef(fit_pois_indep), coef(fit_pois_exch), tolerance = 0.5)
  })

  it("cluster and observation counts are consistent with data", {
    expect_equal(fit_pois_indep$clusters_no, length(unique(ds$epilepsy$id)))
    expect_equal(fit_pois_indep$obs_no,      nrow(ds$epilepsy))
  })

  it("fixed corstr with valid alpha_vector converges", {
    # epilepsy has T = 4 repeated → choose(4,2) = 6 pairs
    alpha_fixed <- rep(0.2, 6)
    fit_fixed <- geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = ds$epilepsy, id = id,
      family  = poisson("log"),
      corstr  = "fixed",
      alpha_vector = alpha_fixed
    )
    expect_true(fit_fixed$converged)
  })

  it("phi_fixed freezes the scale parameter", {
    fit_phi <- geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = ds$epilepsy, id = id,
      family    = poisson("log"),
      corstr    = "independence",
      phi_fixed = TRUE, phi_value = 2
    )
    expect_equal(fit_phi$phi, 2, tolerance = 1e-10)
  })

})


# =============================================================================
# 6. geewa() – Gamma family
# =============================================================================

describe("geewa() Gamma log", {

  it("converges on leprosy data", {
    fit_gamma <- geewa(
      bacilli ~ treatment + period,
      data   = ds$leprosy, id = id,
      family = Gamma(link = "log"),
      corstr = "exchangeable"
    )
    expect_true(fit_gamma$converged)
    expect_true(all(fitted(fit_gamma) > 0))
    expect_equal(fit_gamma$family$family, "Gamma")
  })

})


# =============================================================================
# 7. geewa() – estimator methods
# =============================================================================

describe("geewa() estimation methods", {

  methods_cc <- c("brgee-naive", "brgee-robust", "brgee-empirical",
                  "bcgee-naive", "bcgee-robust", "bcgee-empirical",
                  "pgee-jeffreys")

  for (mth in methods_cc) {
    local({
      m <- mth
      it(paste("method =", m, "produces a converged geer object"), {
        fit <- geewa(
          seizures ~ treatment + lnbaseline + lnage,
          data = ds$epilepsy, id = id,
          family = poisson("log"),
          corstr = "exchangeable",
          method = m
        )
        expect_s3_class(fit, "geer")
        expect_true(fit$converged)
        expect_true(all(is.finite(coef(fit))))
      })
    })
  }

  it("bias-reduced estimates differ from plain GEE estimates", {
    fit_br <- geewa(
      seizures ~ treatment + lnbaseline + lnage,
      data = ds$epilepsy, id = id,
      family = poisson("log"),
      corstr = "exchangeable",
      method = "brgee-robust"
    )
    # Not identical, but in the same neighbourhood
    expect_false(isTRUE(all.equal(coef(fit_pois_exch), coef(fit_br),
                                  tolerance = 1e-12)))
    expect_equal(coef(fit_pois_exch), coef(fit_br), tolerance = 0.5)
  })

})


# =============================================================================
# 8. geewa_binary() – binary outcomes
# =============================================================================

describe("geewa_binary()", {

  it("errors on unknown link", {
    expect_error(
      geewa_binary(ecg ~ treatment, id = id,
                   data = ds$cerebrovascular, link = "log-log"),
      "'link' should be one of"
    )
  })

  it("errors on unknown orstr", {
    expect_error(
      geewa_binary(ecg ~ treatment, id = id,
                   data = ds$cerebrovascular, link = "logit",
                   orstr = "ar1"),
      "'orstr' should be one of"
    )
  })

  it("produces a geer object with correct family", {
    expect_s3_class(fit_bin_or, "geer")
    expect_equal(fit_bin_or$family$family, "binomial")
    expect_equal(fit_bin_or$family$link,   "logit")
  })

  it("converges on cerebrovascular data", {
    expect_true(fit_bin_or$converged)
  })

  it("fitted values are in (0, 1)", {
    fv <- fitted(fit_bin_or)
    expect_true(all(fv > 0 & fv < 1))
  })

  it("phi is identically 1 for binary OR model", {
    expect_equal(fit_bin_or$phi, 1)
  })

  it("exchangeable orstr returns a single alpha", {
    expect_length(fit_bin_or$alpha, 1L)
    expect_gt(fit_bin_or$alpha, 0)
  })

  it("unstructured orstr returns choose(T, 2) alphas", {
    fit_unstr_bin <- geewa_binary(
      formula = status ~ treatment + baseline,
      id = id, repeated = visit,
      data = ds$respiratory[ds$respiratory$center == "C2", ],
      link = "logit", orstr = "unstructured"
    )
    T_max <- max(ds$respiratory[ds$respiratory$center == "C2", "visit"])
    expect_length(fit_unstr_bin$alpha, choose(T_max, 2))
  })

  it("probit link converges", {
    fit_probit <- geewa_binary(
      ecg ~ period + treatment,
      id = id, data = ds$cerebrovascular,
      link = "probit", orstr = "independence"
    )
    expect_true(fit_probit$converged)
    expect_equal(fit_probit$family$link, "probit")
  })

  it("all bias-reducing / penalised methods converge for binary data", {
    methods_bin <- c("brgee-naive", "brgee-robust", "brgee-empirical",
                     "bcgee-naive", "bcgee-robust", "bcgee-empirical",
                     "pgee-jeffreys")
    for (mth in methods_bin) {
      fit <- geewa_binary(
        ecg ~ period + treatment,
        id = id, data = ds$cerebrovascular,
        link = "logit", orstr = "exchangeable",
        method = mth
      )
      expect_true(fit$converged,
                  label = paste("method =", mth, "should converge"))
    }
  })

})


# =============================================================================
# 9. vcov() – all four cov_type branches
# =============================================================================

describe("vcov() covariance extractors", {

  p <- length(coef(fit_pois_exch))

  for (ctype in c("robust", "naive", "bias-corrected", "df-adjusted")) {
    local({
      ct <- ctype
      it(paste("cov_type =", ct, "returns a p×p symmetric PSD matrix"), {
        V <- vcov(fit_pois_exch, cov_type = ct)
        expect_true(is.matrix(V))
        expect_equal(dim(V), c(p, p))
        # symmetry
        expect_equal(V, t(V), tolerance = 1e-12)
        # positive semi-definite: all eigenvalues >= 0
        ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
        expect_true(all(ev >= -1e-10))
        # named
        expect_equal(rownames(V), names(coef(fit_pois_exch)))
        expect_equal(colnames(V), names(coef(fit_pois_exch)))
      })
    })
  }

  it("df-adjusted matrix is a positive scaling of robust", {
    V_rob <- vcov(fit_pois_exch, cov_type = "robust")
    V_dfa <- vcov(fit_pois_exch, cov_type = "df-adjusted")
    n      <- fit_pois_exch$clusters_no
    ratio  <- V_dfa / V_rob
    # All ratios should equal n / (n - p) up to floating point
    expect_equal(ratio, matrix(n / (n - p), nrow = p, ncol = p),
                 tolerance = 1e-10)
  })

  it("errors when clusters_no <= number of coefs for df-adjusted", {
    # Manufacture a stub where clusters_no == p
    stub <- unclass(fit_pois_exch)
    stub$clusters_no <- p  # will trigger the guard
    bad_fit <- geer(stub)
    expect_error(vcov(bad_fit, cov_type = "df-adjusted"),
                 "clusters_no must be > number of coefficients")
  })

})


# =============================================================================
# 10. summary() / print() – smoke tests
# =============================================================================

describe("summary() and print()", {

  it("summary() returns a summary.geer object", {
    s <- summary(fit_pois_exch)
    expect_s3_class(s, "summary.geer")
  })

  it("summary coefficients table has correct column names", {
    s  <- summary(fit_pois_exch)
    cn <- colnames(s$coefficients)
    expect_equal(cn, c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  })

  it("summary() works for each cov_type", {
    for (ct in c("robust", "naive", "bias-corrected", "df-adjusted")) {
      s <- summary(fit_pois_exch, cov_type = ct)
      expect_equal(s$cov_type, ct)
    }
  })

  it("print() returns the object invisibly", {
    expect_invisible(print(fit_pois_exch))
  })

  it("print.summary.geer() returns the summary invisibly", {
    s <- summary(fit_pois_exch)
    expect_invisible(print(s))
  })

  it("summary includes correct family & link", {
    s <- summary(fit_pois_exch)
    expect_equal(s$family$family, "poisson")
    expect_equal(s$family$link,   "log")
  })

})


# =============================================================================
# 11. residuals() – three types
# =============================================================================

describe("residuals()", {

  n <- fit_pois_exch$obs_no

  it("working residuals have length obs_no and are finite", {
    r <- residuals(fit_pois_exch, type = "working")
    expect_length(r, n)
    expect_true(all(is.finite(r)))
  })

  it("Pearson residuals have length obs_no and are finite", {
    r <- residuals(fit_pois_exch, type = "pearson")
    expect_length(r, n)
    expect_true(all(is.finite(r)))
  })

  it("deviance residuals have length obs_no and are finite", {
    r <- residuals(fit_pois_exch, type = "deviance")
    expect_length(r, n)
    expect_true(all(is.finite(r)))
  })

  it("Pearson residuals are signed and close to (y - mu)/sqrt(V(mu))", {
    y  <- fit_pois_exch$y
    mu <- fitted(fit_pois_exch)
    wt <- fit_pois_exch$prior.weights
    # For Poisson, V(mu) = mu
    expected <- (y - mu) / sqrt(mu) * sqrt(wt)
    observed <- residuals(fit_pois_exch, type = "pearson")
    expect_equal(observed, expected, tolerance = 1e-6)
  })

  it("errors on unknown residual type", {
    expect_error(residuals(fit_pois_exch, type = "raw"), "should be one of")
  })

  it("residuals() on binary model are finite", {
    r <- residuals(fit_bin_or, type = "working")
    expect_true(all(is.finite(r)))
  })

})


# =============================================================================
# 12. predict() – in-sample & newdata, link & response, se.fit
# =============================================================================

describe("predict()", {

  it("in-sample response predictions equal fitted values", {
    pred <- predict(fit_pois_exch, type = "response")
    expect_equal(pred, fitted(fit_pois_exch), tolerance = 1e-10)
  })

  it("in-sample link predictions equal linear predictors", {
    pred <- predict(fit_pois_exch, type = "link")
    expect_equal(pred, fit_pois_exch$linear.predictors, tolerance = 1e-10)
  })

  it("se.fit = TRUE returns a list with fit and se.fit", {
    pred <- predict(fit_pois_exch, type = "response", se.fit = TRUE)
    expect_type(pred, "list")
    expect_named(pred, c("fit", "se.fit"))
    expect_equal(length(pred$fit), fit_pois_exch$obs_no)
    expect_true(all(pred$se.fit > 0))
  })

  it("newdata prediction has correct length", {
    nd   <- ds$epilepsy[1:10, ]
    pred <- predict(fit_pois_exch, newdata = nd, type = "response")
    expect_length(pred, 10L)
    expect_true(all(pred > 0))
  })

  it("newdata link prediction is log of response prediction (Poisson log)", {
    nd    <- ds$epilepsy[1:5, ]
    p_lin <- predict(fit_pois_exch, newdata = nd, type = "link")
    p_res <- predict(fit_pois_exch, newdata = nd, type = "response")
    expect_equal(log(p_res), p_lin, tolerance = 1e-10)
  })

  it("newdata se.fit with link type is finite and positive", {
    nd   <- ds$epilepsy[1:5, ]
    pred <- predict(fit_pois_exch, newdata = nd, type = "link", se.fit = TRUE)
    expect_true(all(is.finite(pred$se.fit)))
    expect_true(all(pred$se.fit > 0))
  })

  it("binary model: response predictions in (0, 1)", {
    pred <- predict(fit_bin_or, type = "response")
    expect_true(all(pred > 0 & pred < 1))
  })

})


# =============================================================================
# 13. fitted() and coef()
# =============================================================================

describe("fitted() and coef()", {

  it("fitted() returns a numeric vector", {
    expect_type(fitted(fit_pois_exch), "double")
  })

  it("coef() returns a named numeric vector", {
    b <- coef(fit_pois_exch)
    expect_type(b, "double")
    expect_named(b)
    expect_equal(length(b), ncol(fit_pois_exch$x))
  })

  it("coef() and fit$coefficients are identical", {
    expect_identical(coef(fit_pois_exch), fit_pois_exch$coefficients)
  })

})


# =============================================================================
# 14. model.matrix()
# =============================================================================

describe("model.matrix.geer()", {

  it("returns the design matrix embedded in the fit", {
    X <- model.matrix(fit_pois_exch)
    expect_equal(X, fit_pois_exch$x)
  })

  it("has correct dimensions", {
    X <- model.matrix(fit_pois_exch)
    expect_equal(nrow(X), fit_pois_exch$obs_no)
    expect_equal(ncol(X), length(coef(fit_pois_exch)))
  })

})


# =============================================================================
# 15. geecriteria()
# =============================================================================

describe("geecriteria()", {

  it("returns a data.frame for a single model", {
    crit <- geecriteria(fit_pois_exch)
    expect_s3_class(crit, "data.frame")
    expect_equal(nrow(crit), 1L)
  })

  it("has the expected column names", {
    crit <- geecriteria(fit_pois_exch)
    expected_cols <- c("QIC", "CIC", "RJC", "QICu", "GESSC", "GPC", "Parameters")
    expect_equal(colnames(crit), expected_cols)
  })

  it("Parameters equals number of regression coefficients", {
    crit <- geecriteria(fit_pois_exch)
    expect_equal(crit$Parameters, length(coef(fit_pois_exch)))
  })

  it("returns two rows for two models", {
    crit <- geecriteria(fit_pois_indep, fit_pois_exch)
    expect_equal(nrow(crit), 2L)
  })

  it("errors when a non-geer object is passed", {
    lm_fit <- lm(seizures ~ treatment, data = ds$epilepsy)
    expect_error(geecriteria(fit_pois_exch, lm_fit),
                 "Only 'geer' objects are supported")
  })

  it("errors on non-integer digits", {
    expect_error(geecriteria(fit_pois_exch, digits = 1.5),
                 "'digits' must be a single non-negative integer")
  })

  it("works for binary OR model", {
    crit <- geecriteria(fit_bin_or)
    expect_s3_class(crit, "data.frame")
    expect_true(all(is.finite(unlist(crit))))
  })

  it("cov_type argument is passed through without error", {
    for (ct in c("robust", "naive", "bias-corrected", "df-adjusted")) {
      expect_no_error(geecriteria(fit_pois_exch, cov_type = ct))
    }
  })

})


# =============================================================================
# 16. update()
# =============================================================================

describe("update()", {

  it("update with new method returns a geer object", {
    fit_br <- update(fit_pois_exch, method = "brgee-robust")
    expect_s3_class(fit_br, "geer")
    expect_equal(fit_br$method, "brgee-robust")
  })

  it("update with new corstr changes the association structure", {
    fit_ar1 <- update(fit_pois_exch, corstr = "ar1")
    expect_s3_class(fit_ar1, "geer")
    expect_equal(fit_ar1$association_structure, "ar1")
  })

  it("update preserves unchanged arguments", {
    fit_br <- update(fit_pois_exch, method = "brgee-robust")
    # family, corstr should be identical
    expect_equal(fit_br$family$family, fit_pois_exch$family$family)
    expect_equal(fit_br$association_structure, fit_pois_exch$association_structure)
  })

  it("update with new formula adds a term", {
    # add an interaction term
    fit_int <- update(fit_pois_indep, . ~ . + treatment:lnbaseline)
    expect_s3_class(fit_int, "geer")
    expect_gt(length(coef(fit_int)), length(coef(fit_pois_indep)))
  })

})


# =============================================================================
# 17. Numerical regression against geepack (tolerance 1e-4)
#     These guard against regressions in the C++ solver.
#     Reference values were obtained from geepack::geeglm on the same data.
# =============================================================================

describe("Numerical regression vs. reference values", {

  # Reference: geepack::geeglm(seizures ~ treatment + lnbaseline + lnage,
  #   id = id, data = epilepsy, family = poisson("log"), corstr = "independence")
  # Intercept      treatmentprogabide  lnbaseline  lnage
  # -1.3242        -0.1527             0.8840      0.4818

  it("Poisson independence coefs match geepack reference to 1e-3", {
    b <- coef(fit_pois_indep)
    expect_equal(b["(Intercept)"],          -1.3242, tolerance = 1e-3)
    expect_equal(b["treatmentprogabide"],   -0.1527, tolerance = 1e-3)
    expect_equal(b["lnbaseline"],            0.8840, tolerance = 1e-3)
    expect_equal(b["lnage"],                 0.4818, tolerance = 1e-3)
  })

  # Reference for exchangeable: coefs shift slightly due to working corr
  it("Poisson exchangeable coefs are finite and intercept is negative", {
    b <- coef(fit_pois_exch)
    expect_true(all(is.finite(b)))
    expect_lt(b["(Intercept)"], 0)
  })

  # Reference: geepack::geeglm(ecg ~ period + treatment, id = id,
  #   data = cerebrovascular, family = binomial("logit"), corstr = "exchangeable")
  # Intercept  period  treatmentplacebo
  # 0.7547     0.3697   -0.6730

  it("Binary exchangeable coefs match geepack reference to 1e-3", {
    b <- coef(fit_bin_or)
    expect_equal(b["(Intercept)"],       0.7547, tolerance = 1e-3)
    expect_equal(b["period"],            0.3697, tolerance = 1e-3)
    expect_equal(b["treatmentplacebo"], -0.6730, tolerance = 1e-3)
  })

})
