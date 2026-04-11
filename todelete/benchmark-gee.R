library(gee)
library(MASS)
library(testthat)
devtools::load_all()


strip <- function(x) {
  dimnames(x) <- NULL
  x
}


check_gee_fixed_example <- function(label, ref, fixed, tol = 1e-5) {
  naive_ref <- strip(ref$naive.variance)
  describe(label, {
    it("coefficients agree", {
      expect_equal(drop(ref$coefficients), drop(coef(fixed)), tolerance = tol)
    })
    it("robust vcov agrees", {
      expect_equal(
        strip(ref$robust.variance),
        strip(vcov(fixed, cov_type = "robust")),
        tolerance = tol
      )
    })
    it("naive vcov agrees", {
      expect_equal(
        naive_ref,
        strip(vcov(fixed, cov_type = "naive")),
        tolerance = tol
      )
    })
  })
}


# =============================================================================
# Example 1: warpbreaks
# gaussian | exchangeable
# =============================================================================

data(warpbreaks)
ref1 <- gee(
  breaks ~ tension,
  id = wool,
  data = warpbreaks,
  corstr = "exchangeable"
)
R1 <- ref1$working.correlation
fix1 <- geewa(
  formula = breaks ~ tension,
  id = wool,
  data = warpbreaks,
  family = gaussian("identity"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref1$scale,
  alpha_vector = R1[lower.tri(R1)]
)
check_gee_fixed_example(
  label = "Ex1 warpbreaks gaussian exchangeable fixed from gee",
  ref = ref1,
  fixed = fix1
)

# =============================================================================
# Example 2: warpbreaks
# gaussian | AR-M (Mv = 1)
# =============================================================================

ref2 <- gee(
  breaks ~ tension,
  id = wool,
  data = warpbreaks,
  corstr = "AR-M",
  Mv = 1
)
R2 <- ref2$working.correlation
fix2 <- geewa(
  formula = breaks ~ tension,
  id = wool,
  data = warpbreaks,
  family = gaussian("identity"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref2$scale,
  alpha_vector = R2[lower.tri(R2)]
)
check_gee_fixed_example(
  label = "Ex2 warpbreaks gaussian AR-M(1) fixed from gee",
  ref = ref2,
  fixed = fix2
)

# =============================================================================
# Example 3: OME
# binomial(logit) | exchangeable
# =============================================================================

data(OME)
ref3 <- gee(
  cbind(Correct, Trials - Correct) ~ Loud + Age + OME,
  id = ID,
  data = OME,
  family = binomial,
  corstr = "exchangeable"
)
R3 <- ref3$working.correlation
fix3 <- geewa(
  formula = cbind(Correct, Trials - Correct) ~ Loud + Age + OME,
  id = ID,
  data = OME,
  family = binomial("logit"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref3$scale,
  alpha_vector = R3[lower.tri(R3)]
)
check_gee_fixed_example(
  label = "Ex3 OME binomial/logit exchangeable fixed from gee",
  ref = ref3,
  fixed = fix3
)
