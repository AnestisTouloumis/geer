# =============================================================================
# Benchmark: gee (gee::gee) versus geer::geewa
# =============================================================================

library(gee)
library(testthat)
bench_dir <- if (file.exists("benchmark-helpers.R")) "." else "benchmarks"
devtools::load_all(file.path(bench_dir, ".."))
source(file.path(bench_dir, "benchmark-helpers.R"))

gee_tol <- 1e-10

as_ref_gee <- function(ref) {
  tab <- summary(ref)$coefficients
  list(
    coefficients = ref$coefficients,
    naive = ref$naive.variance,
    robust = ref$robust.variance,
    fitted = ref$fitted.values,
    linear_predictors = ref$linear.predictors,
    residuals_working = ref$residuals,
    summary = list(cov_type = "robust", se = tab[, "Robust S.E."],
                   statistic = tab[, "Robust z"])
  )
}

naive_summary_gee <- function(ref) {
  tab <- summary(ref)$coefficients
  list(summary = list(cov_type = "naive", se = tab[, "Naive S.E."],
                      statistic = tab[, "Naive z"]))
}

fixed_from_gee <- function(args, ref) {
  geewa_fixed(args, ref$working.correlation, ref$scale)
}


# =============================================================================
# Example 1: warpbreaks, gaussian, exchangeable
# =============================================================================

data("warpbreaks", package = "datasets")

ref1 <- quietly(gee(breaks ~ tension, id = wool, data = warpbreaks,
                    corstr = "exchangeable", tol = gee_tol))
fit1 <- fixed_from_gee(list(formula = breaks ~ tension, id = quote(wool),
                            data = warpbreaks,
                            family = gaussian(link = "identity")), ref1)

check_fit("Ex1 warpbreaks gaussian exchangeable", as_ref_gee(ref1), fit1)
check_fit("Ex1 warpbreaks naive summary", naive_summary_gee(ref1), fit1)


# =============================================================================
# Example 2: warpbreaks, gaussian, AR-M with Mv = 1
# =============================================================================

ref2 <- quietly(gee(breaks ~ tension, id = wool, data = warpbreaks,
                    corstr = "AR-M", Mv = 1, tol = gee_tol))
fit2 <- fixed_from_gee(list(formula = breaks ~ tension, id = quote(wool),
                            data = warpbreaks,
                            family = gaussian(link = "identity")), ref2)

check_fit("Ex2 warpbreaks gaussian AR-M(1)", as_ref_gee(ref2), fit2)
check_fit("Ex2 warpbreaks naive summary", naive_summary_gee(ref2), fit2)


# =============================================================================
# Example 3: OME (MASS), binomial(logit), exchangeable
# =============================================================================

data("OME", package = "MASS")

ref3 <- quietly(gee(cbind(Correct, Trials - Correct) ~ Loud + Age + OME,
                    id = ID, data = OME, family = binomial,
                    corstr = "exchangeable", tol = gee_tol))
fit3 <- fixed_from_gee(list(
  formula = cbind(Correct, Trials - Correct) ~ Loud + Age + OME,
  id = quote(ID), data = OME, family = binomial(link = "logit")
), ref3)

ref3_list <- as_ref_gee(ref3)
ref3_list$residuals_working <- NULL
check_fit("Ex3 OME binomial/logit exchangeable", ref3_list, fit3)
check_fit("Ex3 OME naive summary", naive_summary_gee(ref3), fit3)

describe("Ex3 OME residuals", {
  trials <- in_geer_order(OME$Trials, fit3)
  res_geer <- residuals(fit3, type = "working")
  res_gee <- in_geer_order(ref3$residuals, fit3)
  p_geer <- fitted(fit3)
  it("gee y is the success count", {
    expect_agree(OME$Correct, ref3$y)
  })
  it("gee residuals are y - fitted.values", {
    expect_agree(ref3$y - ref3$fitted.values, ref3$residuals)
  })
  it("gee residual = Trials * geer residual + (Trials - 1) * fitted", {
    expect_agree(res_gee, trials * res_geer + (trials - 1) * p_geer)
  })
  it("rows with Trials = 1 agree directly", {
    one <- trials == 1
    expect_true(any(one))
    expect_agree(res_gee[one], res_geer[one])
  })
})
