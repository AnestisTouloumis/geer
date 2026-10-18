# =============================================================================
# Benchmark: BCgee (BCgee) versus geer::geewa(method = "bcgee-robust")
# =============================================================================

library(BCgee)
library(gee)
library(testthat)
bench_dir <- if (file.exists("benchmark-helpers.R")) "." else "benchmarks"
devtools::load_all(file.path(bench_dir, ".."))
source(file.path(bench_dir, "benchmark-helpers.R"))

gee_tol <- 1e-10

as_ref_bcgee <- function(ref) {
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

naive_summary_bcgee <- function(ref) {
  tab <- summary(ref)$coefficients
  list(summary = list(cov_type = "naive", se = tab[, "Naive S.E."],
                      statistic = tab[, "Naive z"]))
}

# geer fits with method = "gee" and "bcgee-robust", alpha and phi fixed at
# the values of the gee() fit.
fixed_pair <- function(args, ref) {
  R <- ref$working.correlation
  list(gee = geewa_fixed(args, R, ref$scale),
       bc = geewa_fixed(c(args, method = "bcgee-robust"), R, ref$scale))
}

check_bias <- function(label, ref_bc, fits) {
  describe(label, {
    it("estimated bias agrees", {
      expect_agree(ref_bc$bias, coef(fits$gee) - coef(fits$bc))
    })
  })
}


# =============================================================================
# Example 1: cerebrovascular deficiency, binomial(logit), exchangeable
# =============================================================================

data("cereb", package = "BCgee")

ref1 <- quietly(gee(y ~ Period + Drug, id = id, data = cereb,
                    family = binomial(logit), corstr = "exchangeable",
                    tol = gee_tol))
ref1_bc <- BCgee(ref1)
fits1 <- fixed_pair(list(formula = y ~ Period + Drug, id = quote(id),
                         data = cereb, family = binomial(link = "logit")),
                    ref1)

check_fit("Ex1 cereb gee() coefficients",
          list(coefficients = ref1$coefficients), fits1$gee)
check_fit("Ex1 cereb BCgee", as_ref_bcgee(ref1_bc), fits1$bc)
check_fit("Ex1 cereb naive summary", naive_summary_bcgee(ref1_bc), fits1$bc)
check_bias("Ex1 cereb", ref1_bc, fits1)


# =============================================================================
# Example 2: epileptic seizures, poisson(log), exchangeable, offset(log(t))
# =============================================================================

data("seizure", package = "BCgee")
seiz.l <- reshape(seizure,
                  varying = list(c("base", "y1", "y2", "y3", "y4")),
                  v.names = "y", times = 0:4, direction = "long")
seiz.l <- seiz.l[order(seiz.l$id, seiz.l$time), ]
seiz.l$t <- ifelse(seiz.l$time == 0, 8, 2)
seiz.l$x <- ifelse(seiz.l$time == 0, 0, 1)

ref2 <- quietly(gee(y ~ offset(log(t)) + x + trt + x:trt, id = id,
                    data = seiz.l, corstr = "exchangeable",
                    family = poisson(log), tol = gee_tol))
ref2_bc <- BCgee(ref2)
fits2 <- fixed_pair(list(formula = y ~ offset(log(t)) + x + trt + x:trt,
                         id = quote(id), data = seiz.l,
                         family = poisson(link = "log")), ref2)


ref2_with_offset <- ref2
ref2_with_offset$fitted.values <- exp(ref2$linear.predictors + log(seiz.l$t))
ref2_bc_offset <- BCgee(ref2_with_offset)

offset2 <- in_geer_order(log(seiz.l$t), fits2$gee)
X2 <- model.matrix(~ x + trt + x:trt, data = seiz.l)

check_fit("Ex2 seizure gee() coefficients",
          list(coefficients = ref2$coefficients), fits2$gee)
describe("Ex2 seizure: offset handling of gee and BCgee", {
  it("gee linear predictors = geer linear predictors - offset", {
    expect_agree(in_geer_order(ref2$linear.predictors, fits2$gee),
                 predict(fits2$gee, type = "link") - offset2)
  })
  it("gee fitted values = exp(geer linear predictors - offset)", {
    expect_agree(in_geer_order(ref2$fitted.values, fits2$gee),
                 exp(predict(fits2$gee, type = "link") - offset2))
  })
  it("BCgee corrected linear predictors = X %*% beta, without the offset", {
    expect_agree(ref2_bc$linear.predictors,
                 drop(X2 %*% ref2_bc$coefficients))
  })
  it("with the offset added, BCgee coefficients agree", {
    expect_agree(ref2_bc_offset$coefficients, coef(fits2$bc))
  })
  it("with the offset added, BCgee bias agrees", {
    expect_agree(ref2_bc_offset$bias, coef(fits2$gee) - coef(fits2$bc))
  })
  it("robust covariance differs: BCgee evaluates it without the offset", {
    expect_known_difference(ref2_bc_offset$robust.variance,
                            vcov(fits2$bc, cov_type = "robust"),
                            "BCgee drops the offset after the correction")
  })
})

ref2b <- quietly(gee(y ~ x + trt + x:trt, id = id, data = seiz.l,
                     corstr = "exchangeable", family = poisson(log),
                     tol = gee_tol))
ref2b_bc <- BCgee(ref2b)
fits2b <- fixed_pair(list(formula = y ~ x + trt + x:trt, id = quote(id),
                          data = seiz.l, family = poisson(link = "log")),
                     ref2b)

check_fit("Ex2 seizure without offset BCgee", as_ref_bcgee(ref2b_bc),
          fits2b$bc)
check_fit("Ex2 seizure without offset naive summary",
          naive_summary_bcgee(ref2b_bc), fits2b$bc)
check_bias("Ex2 seizure without offset", ref2b_bc, fits2b)
