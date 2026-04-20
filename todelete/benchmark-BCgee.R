library("BCgee")
library("gee")
library("geepack")
library("testthat")
devtools::load_all()


tol   <- 1e-5
strip <- function(x) {
  dimnames(x) <- NULL
  x
  }
check_bcgee_example <- function(label, ref_bc, est_bc) {
  describe(label, {
    it("bias-corrected coefficients agree", {
      expect_equal(drop(ref_bc$coefficients), drop(coef(est_bc)),
                   tolerance = tol)
    })
    it("robust vcov agrees", {
      expect_equal(strip(ref_bc$robust.variance), strip(vcov(est_bc, type = "robust")),
                   tolerance = tol)
    })
    it("naive vcov agrees", {
      expect_equal( strip(ref_bc$naive.variance), strip(vcov(est_bc, cov_type = "naive")),
                   tolerance = tol)
    })
  })
}

# =============================================================================
# Example 1: Cerebrovascular deficiency
# binomial(logit) | exchangeable
# Reference: Diggle, Liang, Zeger (1994), p. 153
# =============================================================================

data(cereb, package = "BCgee")
mod1 <- y ~ Period + Drug
ref1_gee <- gee(mod1,
                id = id,
                data = cereb,
                family = binomial("logit"),
                corstr = "exchangeable"
                )
ref1_bc <- BCgee(ref1_gee)
R1 <- ref1_gee$working.correlation
est1_bc <- geewa(mod1,
                 id = id,
                 family = binomial("logit"),
                 data = cereb,
                 method = "bcgee-robust",
                 corstr = "fixed",
                 phi_fixed = TRUE,
                 phi_value = ref1_gee$scale,
                 alpha_vector = R1[lower.tri(R1)])
check_bcgee_example(
  label = "Ex1 cereb binomial/logit exchangeable — BCgee vs geewa bcgee-robust",
  ref_bc = ref1_bc,
  est_bc = est1_bc
)

# =============================================================================
# Example 2: Epileptic seizures
# poisson(log) | exchangeable, with offset
# Reference: Diggle, Liang, Zeger (1994), p. 166
# =============================================================================

data(seizure, package = "geepack")
seiz.l <- reshape(seizure,
                  varying   = list(c("base", "y1", "y2", "y3", "y4")),
                  v.names   = "y",
                  times     = 0:4,
                  direction = "long")
seiz.l <- seiz.l[order(seiz.l$id, seiz.l$time), ]
seiz.l$t <- ifelse(seiz.l$time == 0, 8, 2)
seiz.l$x <- ifelse(seiz.l$time == 0, 0, 1)
#mod2 <- y ~ offset(log(t)) + x * trt
# offset term does not work possible bug!? in BCgee
mod2 <- y ~  x + trt + x:trt
ref2_gee <- gee(mod2,
                id = id,
                data = seiz.l,
                corstr = "exchangeable",
                family = poisson("log"))
ref2_bc <- BCgee(ref2_gee)
R2 <- ref2_gee$working.correlation
est2_bc <- geewa(mod2,
                  id = id,
                  family = poisson("log"),
                  data = seiz.l,
                  method = "bcgee-robust",
                  corstr = "fixed",
                  phi_fixed = TRUE,
                  phi_value = ref2_gee$scale,
                  alpha_vector = R2[lower.tri(R2)])

check_bcgee_example(
  label = "Ex2 seizures poisson/log exchangeable — BCgee vs geewa bcgee-robust",
  ref_bc = ref2_bc,
  est_bc = est2_bc
)
