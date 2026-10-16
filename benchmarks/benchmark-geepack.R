# =============================================================================
# Benchmark: geepack (geeglm, geese) versus geer::geewa
# =============================================================================

library(geepack)
library(testthat)
bench_dir <- if (file.exists("benchmark-helpers.R")) "." else "benchmarks"
devtools::load_all(file.path(bench_dir, ".."))
source(file.path(bench_dir, "benchmark-helpers.R"))
source(file.path(bench_dir, "extract-correlation.R"))

gp_control <- geese.control(epsilon = 1e-10, maxit = 100)

run_id <- function(id) cumsum(c(TRUE, id[-1L] != id[-length(id)]))

# Working correlation matrix from the alpha of a geese object, or of the
# geese component of a geeglm object. "toeplitz" has one parameter per lag.
geese_corr <- function(alpha, structure, size) {
  alpha <- unname(alpha)
  lag <- abs(outer(seq_len(size), seq_len(size), "-"))
  R <- diag(size)
  if (structure == "exchangeable") R[lag > 0] <- alpha[1L]
  if (structure == "ar1") R <- alpha[1L]^lag
  if (structure == "toeplitz") {
    for (k in seq_along(alpha)) R[lag == k] <- alpha[k]
  }
  if (structure == "unstructured") {
    R[lower.tri(R)] <- alpha
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
  }
  R
}

tidy_confint <- function(ref) {
  td <- tidy(ref, conf.int = TRUE)
  list(cov_type = "robust", level = 0.95,
       ci = cbind(td$conf.low, td$conf.high))
}

as_ref_geeglm <- function(ref) {
  tab <- coef(summary(ref))
  list(
    coefficients = coef(ref),
    naive = ref$geese$vbeta.naiv,
    robust = vcov(ref),
    fitted = fitted(ref),
    linear_predictors = ref$linear.predictors,
    residuals_working = residuals(ref, type = "response"),
    summary = list(cov_type = "robust", se = tab[, "Std.err"],
                   statistic = tab[, "Wald"], statistic_type = "wald",
                   p_value = tab[, "Pr(>|W|)"]),
    confint = tidy_confint(ref)
  )
}

as_ref_geese <- function(ref) {
  list(coefficients = ref$beta, naive = ref$vbeta.naiv, robust = ref$vbeta)
}

check_pearson_geepack <- function(label, ref, fit) {
  describe(paste(label, "Pearson residuals"), {
    it("geepack residual = sqrt(phi) * geer residual", {
      expect_agree(in_geer_order(residuals(ref, type = "pearson"), fit),
                   sqrt(fit$phi) * residuals(fit, type = "pearson"))
    })
  })
}

# fit_ind: geewa() fit of the working-independence model at the scale of the
# geepack independence refit made by QIC().
check_qic_geepack <- function(label, ref, fit, fit_ind) {
  q <- geepack::QIC(ref)
  p <- length(coef(fit))
  K <- length(unique(ref$id))
  quasi <- fit$phi * (p - geer_criteria(fit, "QICu")[["QICu"]] / 2)
  cic <- sum(diag(solve(vcov(fit_ind, cov_type = "naive"),
                        vcov(fit, cov_type = "robust"))))
  k <- p + length(ref$geese$alpha)
  qic <- -2 * quasi + 2 * cic
  describe(paste(label, "QIC()"), {
    it("Quasi Lik = phi * geer quasi-likelihood", {
      expect_agree(q[["Quasi Lik"]], quasi)
    })
    it("CIC = tr(Omega_I V)", {
      expect_agree(q[["CIC"]], cic)
    })
    it("QIC = -2 Quasi Lik + 2 CIC", {
      expect_agree(q[["QIC"]], qic)
    })
    it("QICu = -2 Quasi Lik + 2 p", {
      expect_agree(q[["QICu"]], -2 * quasi + 2 * p)
    })
    it("QICC = QIC + 2 k (k + 1) / (K - k - 1), K = distinct id values", {
      expect_agree(q[["QICC"]], qic + 2 * k * (k + 1) / (K - k - 1))
    })
    if (K != fit$clusters_no) {
      it("geepack K differs from the number of clusters", {
        expect_known_difference(
          K, fit$clusters_no,
          "QICC of geepack counts distinct id values, not clusters"
        )
      })
    }
  })
}

# geewa() counterpart of the working-independence refit made by QIC().
independence_fit <- function(ref, args, size) {
  ref_ind <- update(ref, corstr = "independence", zcor = NULL)
  geewa_fixed(args, diag(size), ref_ind$geese$gamma)
}

# Fixed geewa() fit from a geeglm object.
fixed_from_geeglm <- function(args, ref, size) {
  geewa_fixed(args, extract_corr_matrix(ref, size), ref$geese$gamma)
}


# =============================================================================
# ?QIC and ?ohio: ohio, binomial(logit), exchangeable and ar1, scale.fix
# =============================================================================

data("ohio", package = "geepack")
args_ohio <- list(formula = resp ~ age + smoke + age:smoke, id = quote(id),
                  data = ohio, family = binomial(link = "logit"))

ref_ohio_ex <- geeglm(resp ~ age + smoke + age:smoke, id = id, data = ohio,
                      family = binomial, corstr = "exch", scale.fix = TRUE,
                      control = gp_control)
ref_ohio_ar <- geeglm(resp ~ age + smoke + age:smoke, id = id, data = ohio,
                      family = binomial, corstr = "ar1", scale.fix = TRUE,
                      control = gp_control)
fit_ohio_ex <- fixed_from_geeglm(args_ohio, ref_ohio_ex, 4)
fit_ohio_ar <- fixed_from_geeglm(args_ohio, ref_ohio_ar, 4)
fit_ohio_ind <- independence_fit(ref_ohio_ex, args_ohio, 4)

check_fit("ohio exchangeable, scale.fix", as_ref_geeglm(ref_ohio_ex),
          fit_ohio_ex)
check_pearson_geepack("ohio exchangeable", ref_ohio_ex, fit_ohio_ex)
check_qic_geepack("ohio exchangeable", ref_ohio_ex, fit_ohio_ex,
                  fit_ohio_ind)
check_fit("ohio ar1, scale.fix", as_ref_geeglm(ref_ohio_ar), fit_ohio_ar)
check_pearson_geepack("ohio ar1", ref_ohio_ar, fit_ohio_ar)
check_qic_geepack("ohio ar1", ref_ohio_ar, fit_ohio_ar, fit_ohio_ind)


# =============================================================================
# ?compCoef: simulated data, gaussian, unstructured (geese)
# =============================================================================

set.seed(1)
gendat_cc <- function(ncl, clsz) {
  id <- rep(1:ncl, each = clsz)
  visit <- rep(1:clsz, ncl)
  n <- ncl * clsz
  x1 <- rbinom(n, 1, 0.5)
  x2 <- runif(n, 0, 1)
  rho <- 2 / 3
  rhomat <- rho^outer(1:4, 1:4, function(x, y) abs(x - y))
  chol.u <- chol(rhomat)
  noise <- as.vector(sapply(1:ncl, function(x) chol.u %*% rnorm(clsz)))
  y <- 1 + 3 * x1 - 2 * x2 + noise
  data.frame(y, id, visit, x1, x2)
}
simdat_cc <- gendat_cc(100, 4)

ref_cc0 <- geese(y ~ x1, id = id, data = simdat_cc, corstr = "un",
                 control = gp_control)
ref_cc1 <- geese(y ~ x1 + x2, id = id, data = simdat_cc, corstr = "un",
                 control = gp_control)
args_cc <- list(id = quote(id), data = simdat_cc)
fit_cc0 <- geewa_fixed(c(list(formula = y ~ x1), args_cc),
                       geese_corr(ref_cc0$alpha, "unstructured", 4),
                       ref_cc0$gamma)
fit_cc1 <- geewa_fixed(c(list(formula = y ~ x1 + x2), args_cc),
                       geese_corr(ref_cc1$alpha, "unstructured", 4),
                       ref_cc1$gamma)

check_fit("compCoef fit0 gaussian unstructured", as_ref_geese(ref_cc0),
          fit_cc0)
check_fit("compCoef fit1 gaussian unstructured", as_ref_geese(ref_cc1),
          fit_cc1)
not_comparable("compCoef(fit0, fit1)",
               "geer has no function comparing coefficients across models")


# =============================================================================
# ?fixed2Zcor: simulated data, gaussian, fixed correlation, permuted waves
# =============================================================================

set.seed(1)
timeorder <- rep(1:5, 6)
tvar <- timeorder + rnorm(length(timeorder))
idvar <- rep(1:6, each = 5)
uuu <- rep(rnorm(6), each = 5)
yvar <- 1 + 2 * tvar + uuu + rnorm(length(tvar))
simdat_fz <- data.frame(idvar, timeorder, tvar, yvar)
simdat_perm <- simdat_fz[sample(nrow(simdat_fz)), ]
simdat_perm <- simdat_perm[order(simdat_perm$idvar), ]
cor_fixed <- matrix(c(1, 0.5, 0.25, 0.125, 0.125,
                      0.5, 1, 0.25, 0.125, 0.125,
                      0.25, 0.25, 1, 0.5, 0.125,
                      0.125, 0.125, 0.5, 1, 0.125,
                      0.125, 0.125, 0.125, 0.125, 1), nrow = 5, ncol = 5)
zcor_fz <- fixed2Zcor(cor_fixed, id = simdat_perm$idvar,
                      waves = simdat_perm$timeorder)
ref_fz <- geeglm(yvar ~ tvar, id = idvar, data = simdat_perm,
                 corstr = "fixed", zcor = zcor_fz, control = gp_control)
fit_fz <- geewa_fixed(list(formula = yvar ~ tvar, id = quote(idvar),
                           repeated = quote(timeorder), data = simdat_perm),
                      cor_fixed, ref_fz$geese$gamma)

check_fit("fixed2Zcor gaussian fixed correlation with waves",
          as_ref_geeglm(ref_fz), fit_fz)
check_pearson_geepack("fixed2Zcor", ref_fz, fit_fz)


# =============================================================================
# ?geeglm: dietox, poisson(identity), ar1, and anova()
# =============================================================================

data("dietox", package = "geepack")
dietox$Cu <- as.factor(dietox$Cu)
args_dietox <- list(data = dietox, id = quote(Pig),
                    family = poisson("identity"))
geeglm_dietox <- function(formula) {
  muffle_noninteger(geeglm(formula, data = dietox, id = Pig,
                           family = poisson("identity"), corstr = "ar1",
                           control = gp_control))
}
fixed_dietox <- function(formula, ref) {
  muffle_noninteger(fixed_from_geeglm(c(list(formula = formula), args_dietox),
                                      ref, 12))
}

mf <- formula(Weight ~ Cu * (Time + I(Time^2) + I(Time^3)))
ref_gee1 <- geeglm_dietox(mf)
fit_gee1 <- fixed_dietox(mf, ref_gee1)
check_fit("geeglm dietox poisson/identity ar1", as_ref_geeglm(ref_gee1),
          fit_gee1)
check_pearson_geepack("geeglm dietox", ref_gee1, fit_gee1)

mf2 <- formula(Weight ~ Cu * Time + I(Time^2) + I(Time^3))
ref_gee2 <- geeglm_dietox(mf2)
fit_gee2 <- fixed_dietox(mf2, ref_gee2)
check_fit("geeglm dietox gee2", as_ref_geeglm(ref_gee2), fit_gee2)

anova_gp <- muffle_noninteger(anova(ref_gee2))
term_labels <- attr(terms(mf2), "term.labels")
seq_formulas <- c(list(Weight ~ 1),
                  lapply(seq_along(term_labels), function(i) {
                    reformulate(term_labels[seq_len(i)], response = "Weight")
                  }))
anova_rel <- sapply(seq_along(term_labels), function(i) {
  ref_i <- geeglm_dietox(seq_formulas[[i + 1L]])
  small <- fixed_dietox(seq_formulas[[i]], ref_i)
  large <- fixed_dietox(seq_formulas[[i + 1L]], ref_i)
  a <- muffle_noninteger(anova(small, large, test = "wald",
                               cov_type = "robust"))
  c(Df = a$Df[2L], Chi = a$Chi[2L], error = ref_i$geese$error)
})
converged <- anova_rel["error", ] == 0
describe("geeglm dietox gee2 anova()", {
  it("degrees of freedom agree", {
    expect_agree(anova_gp$Df, anova_rel["Df", ])
  })
  it("Wald statistics equal geer anova() of the neighbouring models", {
    expect_agree(anova_gp$X2[converged], anova_rel["Chi", converged])
  })
  it("rows whose geepack refit did not converge differ", {
    expect_known_difference(anova_gp$X2[!converged],
                            anova_rel["Chi", !converged],
                            "the geepack refit stopped with error code 1")
  })
})


# =============================================================================
# ?geese and ?seizure: seizure, poisson(log), exchangeable and fixed
# =============================================================================

data("seizure", package = "geepack")
seiz.l <- reshape(seizure,
                  varying = list(c("base", "y1", "y2", "y3", "y4")),
                  v.names = "y", times = 0:4, direction = "long")
seiz.l <- seiz.l[order(seiz.l$id, seiz.l$time), ]
seiz.l$t <- ifelse(seiz.l$time == 0, 8, 2)
seiz.l$x <- ifelse(seiz.l$time == 0, 0, 1)
args_seiz <- list(formula = y ~ offset(log(t)) + x + trt + x:trt,
                  id = quote(id), data = seiz.l, family = poisson)

ref_m1 <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
                data = seiz.l, corstr = "exch", family = poisson,
                control = gp_control)
fit_m1 <- geewa_fixed(args_seiz, geese_corr(ref_m1$alpha, "exchangeable", 5),
                      ref_m1$gamma)
check_fit("geese seizure m1 poisson exchangeable with offset",
          as_ref_geese(ref_m1), fit_m1)

ref_m2 <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
                data = seiz.l, subset = id != 49, corstr = "exch",
                family = poisson, control = gp_control)
fit_m2 <- geewa_fixed(c(args_seiz, subset = quote(id != 49)),
                      geese_corr(ref_m2$alpha, "exchangeable", 5),
                      ref_m2$gamma)
check_fit("geese seizure m2 with subset = id != 49", as_ref_geese(ref_m2),
          fit_m2)

zcor_m3 <- rep(cor_fixed[lower.tri(cor_fixed)], 59)
ref_m3 <- geese(y ~ offset(log(t)) + x + trt + x:trt, id = id,
                data = seiz.l, family = poisson, corstr = "fixed",
                zcor = zcor_m3, control = gp_control)
fit_m3 <- geewa_fixed(args_seiz, cor_fixed, ref_m3$gamma)
check_fit("geese seizure m3 fixed correlation", as_ref_geese(ref_m3), fit_m3)

ref_geese_ex <- geese(resp ~ age + smoke + age:smoke, id = id, data = ohio,
                      family = binomial, corstr = "exch", scale.fix = TRUE,
                      control = gp_control)
ref_geese_ar <- geese(resp ~ age + smoke + age:smoke, id = id, data = ohio,
                      family = binomial, corstr = "ar1", scale.fix = TRUE,
                      control = gp_control)
check_fit("geese ohio exchangeable, scale.fix", as_ref_geese(ref_geese_ex),
          geewa_fixed(args_ohio,
                      geese_corr(ref_geese_ex$alpha, "exchangeable", 4),
                      ref_geese_ex$gamma))
check_fit("geese ohio ar1, scale.fix", as_ref_geese(ref_geese_ar),
          geewa_fixed(args_ohio, geese_corr(ref_geese_ar$alpha, "ar1", 4),
                      ref_geese_ar$gamma))

not_comparable("geese simulated data with sformula = ~ x1",
               "geer has no model for the scale parameter (sformula)")
not_comparable("geese userdefined zcor with sformula and jackknife",
               "scale model (sformula) and jack/j1s/fij variances")
not_comparable("geese userdefined zcor with tied lag-1 correlations",
               "geer has no user-defined correlation design (zcor)")
not_comparable("geese non-constant cluster sizes with sformula",
               "geer has no model for the scale parameter (sformula)")
not_comparable("seizure m3 and m4 (Thall and Vail) with sformula",
               "scale model (sformula), cor.link and sca.link")


# =============================================================================
# ?genZcor: Toeplitz structure through a userdefined zcor
# =============================================================================

gendat_tz <- function() {
  id <- gl(5, 4, 20)
  visit <- rep(1:4, 5)
  y <- rnorm(id)
  data.frame(y, id, visit)[c(-2, -9), ]
}
set.seed(88)
dat_tz <- gendat_tz()
zcor_tz <- genZcor(clusz = table(dat_tz$id), waves = dat_tz$visit,
                   corstrv = 4)
zcor_toep <- matrix(NA, nrow(zcor_tz), 3)
zcor_toep[, 1] <- apply(zcor_tz[, c(1, 4, 6)], 1, sum)
zcor_toep[, 2] <- apply(zcor_tz[, c(2, 5)], 1, sum)
zcor_toep[, 3] <- zcor_tz[, 3]
ref_zfit1 <- geese(y ~ 1, id = id, data = dat_tz, corstr = "userdefined",
                   zcor = zcor_toep, control = gp_control)
ref_zfit2 <- geeglm(y ~ 1, id = id, data = dat_tz, corstr = "userdefined",
                    zcor = zcor_toep, control = gp_control)
args_tz <- list(formula = y ~ 1, id = quote(id), repeated = quote(visit),
                data = dat_tz)
fit_zfit1 <- geewa_fixed(args_tz, geese_corr(ref_zfit1$alpha, "toeplitz", 4),
                         ref_zfit1$gamma)
fit_zfit2 <- geewa_fixed(args_tz,
                         geese_corr(ref_zfit2$geese$alpha, "toeplitz", 4),
                         ref_zfit2$geese$gamma)

check_fit("genZcor geese Toeplitz", as_ref_geese(ref_zfit1), fit_zfit1)
check_fit("genZcor geeglm Toeplitz", as_ref_geeglm(ref_zfit2), fit_zfit2)
check_pearson_geepack("genZcor geeglm", ref_zfit2, fit_zfit2)


# =============================================================================
# ?ordgee: ordinal GEE with an odds-ratio association
#
# For a binary response, ordgee() is GEE for the marginal logit of
# P(y = 1) with the association modelled by log odds ratios, which is
# geewa_binary(). Its alpha holds the log odds ratios, so geewa_binary() is
# fitted with orstr = "fixed" and alpha_vector = exp(alpha) (all 1 for
# "independence"). The koch and respdis examples have three response
# categories, which geer does not model.
# =============================================================================

data("ohio", package = "geepack")
ohio_ord <- ohio
ohio_ord$resp <- ordered(as.factor(ohio_ord$resp))

as_ref_ordgee <- function(ref) {
  tab <- summary(ref)$mean
  list(coefficients = ref$beta, naive = ref$vbeta.naiv, robust = ref$vbeta,
       summary = list(cov_type = "robust", se = tab[, "san.se"],
                      statistic = tab[, "wald"], statistic_type = "wald",
                      p_value = tab[, "p"]))
}

fixed_from_ordgee <- function(ref) {
  log_or <- if (length(ref$alpha) > 0) ref$alpha else 0
  geewa_binary(formula = resp ~ age + smoke + age:smoke, id = id,
               data = ohio, orstr = "fixed",
               alpha_vector = rep(exp(unname(log_or)), length.out = 6),
               control = bench_control)
}

ref_ord <- ordgee(resp ~ age + smoke + age:smoke, id = id, data = ohio_ord,
                  control = gp_control)
check_fit("ordgee ohio binary, independence", as_ref_ordgee(ref_ord),
          fixed_from_ordgee(ref_ord))

# The same model with the other association structures of ordgee() (not
# part of the example).
for (cs in c("exchangeable", "unstructured")) {
  ref_ord_cs <- ordgee(resp ~ age + smoke + age:smoke, id = id,
                       data = ohio_ord, corstr = cs, control = gp_control)
  check_fit(paste("ordgee ohio binary,", cs), as_ref_ordgee(ref_ord_cs),
            fixed_from_ordgee(ref_ord_cs))
}

not_comparable("koch and respdis examples (ordgee)",
               "three ordered response categories; geer models binary data")


# =============================================================================
# ?muscatine: binomial(logit), independence, waves = occasion
# =============================================================================

data("muscatine", package = "geepack")
muscatine$cage <- muscatine$age - 12
muscatine$cage2 <- muscatine$cage^2
f1 <- numobese ~ gender
f2 <- numobese ~ gender + cage + cage2 + gender:cage + gender:cage2

for (f in list(f1, f2)) {
  ref_musc <- geeglm(formula = f, id = id, waves = occasion,
                     data = muscatine, family = binomial(),
                     corstr = "independence", control = gp_control)
  args_musc <- list(formula = f, id = quote(id), repeated = quote(occasion),
                    data = muscatine, family = binomial())
  fit_musc <- fixed_from_geeglm(args_musc, ref_musc, 3)
  label <- paste("muscatine", deparse(f[[3]]))
  check_fit(label, as_ref_geeglm(ref_musc), fit_musc)
  check_pearson_geepack(label, ref_musc, fit_musc)
  check_qic_geepack(label, ref_musc, fit_musc, fit_musc)
}


# =============================================================================
# ?respiratory: binomial(logit), four working structures
# =============================================================================

data("respiratory", package = "geepack")
respiratory$center <- factor(respiratory$center)
respiratory$run_id <- run_id(respiratory$id)
mf_resp <- outcome ~ center + treat + age + baseline
args_resp <- list(formula = mf_resp, data = respiratory, id = quote(run_id),
                  family = binomial())

ref_resp <- list()
for (cs in c("independence", "exchangeable", "unstructured", "ar1")) {
  ref_resp[[cs]] <- geeglm(mf_resp, data = respiratory, id = id,
                           family = binomial(), corstr = cs,
                           control = gp_control)
}
fit_resp_ind <- fixed_from_geeglm(args_resp, ref_resp[["independence"]], 4)
for (cs in names(ref_resp)) {
  fit_resp <- fixed_from_geeglm(args_resp, ref_resp[[cs]], 4)
  label <- paste("respiratory", cs)
  check_fit(label, as_ref_geeglm(ref_resp[[cs]]), fit_resp)
  check_pearson_geepack(label, ref_resp[[cs]], fit_resp)
  check_qic_geepack(label, ref_resp[[cs]], fit_resp, fit_resp_ind)
}

ref_glm <- glm(mf_resp, data = respiratory, family = binomial())
describe("respiratory glm()", {
  it("coefficients equal those of the independence fit", {
    expect_agree(coef(ref_glm), coef(fit_resp_ind))
  })
})


# =============================================================================
# ?relRisk: relative risk regression by the COPY method
# =============================================================================

not_comparable("relRisk respiratory ar1",
               "relRisk fits a log-binomial model by the COPY method")


# =============================================================================
# ?spruce: gaussian, ar1 (geese), 1988 and 1989 seasons
# =============================================================================

data("spruce", package = "geepack")
spruce$contr <- ifelse(spruce$ozone == "enriched", 0, 1)
sitka88 <- spruce[spruce$wave <= 5, ]
sitka89 <- spruce[spruce$wave > 5, ]

ref_88 <- geese(logsize ~ as.factor(wave) + contr + I(time / 100 * contr) - 1,
                id = id, data = sitka88, corstr = "ar1", control = gp_control)
ref_89 <- geese(logsize ~ as.factor(wave) + contr - 1,
                id = id, data = sitka89, corstr = "ar1", control = gp_control)
fit_88 <- geewa_fixed(
  list(formula = logsize ~ as.factor(wave) + contr +
         I(time / 100 * contr) - 1,
       id = quote(id), data = sitka88),
  geese_corr(ref_88$alpha, "ar1", 5), ref_88$gamma
)
fit_89 <- geewa_fixed(
  list(formula = logsize ~ as.factor(wave) + contr - 1,
       id = quote(id), data = sitka89),
  geese_corr(ref_89$alpha, "ar1", 8), ref_89$gamma
)

check_fit("spruce 1988 gaussian ar1", as_ref_geese(ref_88), fit_88)
check_fit("spruce 1989 gaussian ar1", as_ref_geese(ref_89), fit_89)
