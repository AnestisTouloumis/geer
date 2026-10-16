# =============================================================================
# Benchmark: geessbin (geessbin) versus geer::geewa
# =============================================================================

library(geessbin)
library(testthat)
bench_dir <- if (file.exists("benchmark-helpers.R")) "." else "benchmarks"
devtools::load_all(file.path(bench_dir, ".."))
source(file.path(bench_dir, "benchmark-helpers.R"))

beta_map <- c(GEE = "gee", BCGEE = "bcgee-robust", PGEE = "pgee-jeffreys")
se_map <- c(SA = "robust", MK = "df-adjusted", MB = "bias-corrected")

# geessbin's Morel et al. covariance, from the geer naive and robust matrices.
morel_kappa_in_trace <- function(fit) {
  K <- fit$clusters_no
  N <- fit$obs_no
  p <- length(coef(fit))
  naive <- vcov(fit, cov_type = "naive")
  robust <- vcov(fit, cov_type = "robust")
  meat <- solve(naive, t(solve(naive, robust)))
  kappa <- (N - 1) / (N - p) * K / (K - 1)
  lambda <- min(0.5, p / (K - p))
  xi <- max(1, sum(diag(naive %*% (kappa * meat))) / p)
  kappa * robust + lambda * xi * naive
}

data("wheeze", package = "geessbin")

# geessbin() with uncentered cluster scores in the MB covariance.
geessbin_uncentered <- geessbin::geessbin
body(geessbin_uncentered) <- parse(text = sub(
  "cov(d)", "crossprod(d) / (K - 1)",
  paste(deparse(body(geessbin::geessbin)), collapse = "\n"), fixed = TRUE
))[[1L]]
environment(geessbin_uncentered) <- asNamespace("geessbin")

fit_geessbin <- function(corstr, beta_method, se_method,
                         fun = geessbin::geessbin) {
  fun(formula = Wheeze ~ City + factor(Age), data = wheeze, id = ID,
      corstr = corstr, repeated = Age, beta.method = beta_method,
      SE.method = se_method, tol = 1e-10, maxitr = 200)
}

fixed_from_geessbin <- function(ref) {
  geewa_fixed(list(formula = Wheeze ~ City + factor(Age), data = wheeze,
                   id = quote(ID), repeated = quote(Age),
                   family = binomial(link = "logit"),
                   method = beta_map[[ref$beta.method]]),
              ref$wcorr, ref$scale)
}

check_geessbin <- function(label, ref, fit) {
  cov_type <- se_map[[ref$SE.method]]
  s <- summary(ref)
  r <- list(
    coefficients = ref$coefficients,
    fitted = ref$fitted.values,
    linear_predictors = ref$linear.predictors,
    residuals_working = ref$residuals,
    ordered = TRUE
  )
  if (cov_type != "bias-corrected") {
    r[[sub("-", "_", cov_type)]] <- ref$covb
    r$summary <- list(cov_type = cov_type,
                      se = s$coefficients[, "Std.err"],
                      statistic = s$coefficients[, "Z"],
                      p_value = s$coefficients[, "Pr(>|Z|)"])
  }
  check_fit(label, r, fit)
  describe(paste(label, "relationships"), {
    if (cov_type != "bias-corrected") {
      it("odds-ratio limits = exp(confint())", {
        keep <- names(coef(fit)) != "(Intercept)"
        ci <- confint(fit, level = ref$conf.level, cov_type = cov_type)
        expect_agree(s$OR[, c("Lower Limit", "Upper Limit")],
                     exp(ci[keep, ]))
      })
    } else if (identical(ref$beta.method, "GEE")) {
      it("MB = Morel et al. with kappa inside xi", {
        expect_agree(ref$covb, morel_kappa_in_trace(fit))
      })
    } else {
      it("MB differs: geessbin centers the cluster scores", {
        expect_known_difference(
          ref$covb, morel_kappa_in_trace(fit),
          "geessbin centers the cluster scores, which geer does not return"
        )
      })
      it("uncentered MB = Morel et al. with kappa inside xi", {
        ref_u <- fit_geessbin(ref$corstr, ref$beta.method, "MB",
                              fun = geessbin_uncentered)
        expect_agree(ref_u$covb, morel_kappa_in_trace(fit))
      })
    }
  })
}


# =============================================================================
# Example: wheeze, PGEE, ar1, Morel et al. standard errors
# =============================================================================

ref_ex <- fit_geessbin("ar1", "PGEE", "MB")
check_geessbin("Example wheeze ar1 PGEE MB", ref_ex,
               fixed_from_geessbin(ref_ex))


# =============================================================================
# All combinations of working structure, estimator and SE estimator
# =============================================================================

for (corstr in c("independence", "exchangeable", "ar1", "unstructured")) {
  for (beta_method in names(beta_map)) {
    for (se_method in names(se_map)) {
      ref <- fit_geessbin(corstr, beta_method, se_method)
      check_geessbin(paste("wheeze", corstr, beta_method, se_method), ref,
                     fixed_from_geessbin(ref))
    }
  }
}

not_comparable("SE.method KC, MD, FG, PA, GS, WL, WB, FW and FZ",
               "geer has no counterpart of these covariance estimators")
