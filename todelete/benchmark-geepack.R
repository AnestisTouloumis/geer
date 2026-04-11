library(geepack)
library(testthat)
devtools::load_all()

tol <- 1e-5

strip <- function(x) {
  dimnames(x) <- NULL
  x
}

check_geepack_fixed_example <- function(label, ref, fixed, max_size, tol = 1e-5) {
  R <- extract_corr_matrix(ref, max_size)
  naive_ref <- strip(ref$geese$vbeta.naiv)
  describe(label, {
    it("coefficients agree", {
      expect_equal(drop(coef(ref)), drop(coef(fixed)), tolerance = tol)
    })
    it("robust vcov agrees", {
      expect_equal(
        strip(vcov(ref)),
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
# Example 1: dietox - poisson(identity), ar1
# Weight ~ Cu * (Time + I(Time^2) + I(Time^3))
# =============================================================================

data(dietox)
dietox$Cu <- as.factor(dietox$Cu)
mf1 <- formula(Weight ~ Cu * (Time + I(Time^2) + I(Time^3)))

ref1 <- geeglm(
  mf1,
  data = dietox,
  id = Pig,
  family = poisson("identity"),
  corstr = "ar1"
)
R1 <- extract_corr_matrix(ref1, 12)
fix1 <- geewa(
  formula = mf1,
  data = dietox,
  id = Pig,
  family = poisson("identity"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref1$geese$gamma,
  alpha_vector = R1[lower.tri(R1)]
)

check_geepack_fixed_example(
  label = "Ex1 dietox poisson/identity ar1 fixed from geepack",
  ref = ref1,
  fixed = fix1,
  max_size = 12
)

# =============================================================================
# Example 2: ChickWeight - gaussian(identity), independence
# Data must be sorted by Chick for contiguous clusters.
# =============================================================================

chick1 <- ChickWeight[order(ChickWeight$Chick), ]

ref2 <- geeglm(
  weight ~ Time,
  id = Chick,
  data = chick1
)
R2 <- extract_corr_matrix(ref2, 12)
fix2 <- geewa(
  formula = weight ~ Time,
  id = Chick,
  data = chick1,
  family = gaussian("identity"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref2$geese$gamma,
  alpha_vector = R2[lower.tri(R2)]
)

check_geepack_fixed_example(
  label = "Ex2 ChickWeight gaussian/identity independence fixed from geepack",
  ref = ref2,
  fixed = fix2,
  max_size = 12
)

# =============================================================================
# Example 3: ohio - binomial, exchangeable, scale.fix = TRUE (phi = 1)
# =============================================================================

data(ohio)

ref3 <- geeglm(
  resp ~ age + smoke + age:smoke,
  id = id,
  data = ohio,
  family = binomial,
  corstr = "exch",
  scale.fix = TRUE
)
R3 <- extract_corr_matrix(ref3, 4)
fix3 <- geewa(
  formula = resp ~ age + smoke + age:smoke,
  id = id,
  data = ohio,
  family = binomial("logit"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref3$geese$gamma,
  alpha_vector = R3[lower.tri(R3)]
)

check_geepack_fixed_example(
  label = "Ex3 ohio binomial/logit exchangeable fixed from geepack",
  ref = ref3,
  fixed = fix3,
  max_size = 4
)

# =============================================================================
# Example 4: ohio - binomial, ar1, scale.fix = TRUE (phi = 1)
# =============================================================================

ref4 <- geeglm(
  resp ~ age + smoke + age:smoke,
  id = id,
  data = ohio,
  family = binomial,
  corstr = "ar1",
  scale.fix = TRUE
)
R4 <- extract_corr_matrix(ref4, 4)
fix4 <- geewa(
  formula = resp ~ age + smoke + age:smoke,
  id = id,
  data = ohio,
  family = binomial("logit"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref4$geese$gamma,
  alpha_vector = R4[lower.tri(R4)]
)

check_geepack_fixed_example(
  label = "Ex4 ohio binomial/logit ar1 fixed from geepack",
  ref = ref4,
  fixed = fix4,
  max_size = 4
)

# =============================================================================
# Example 5: muscatine - binomial, independence (simple model)
# numobese ~ gender
# =============================================================================

data(muscatine)
muscatine$cage  <- muscatine$age - 12
muscatine$cage2 <- muscatine$cage^2

ref5 <- geeglm(
  numobese ~ gender,
  id = id,
  waves = occasion,
  data = muscatine,
  family = binomial(),
  corstr = "independence"
)
R5 <- extract_corr_matrix(ref5, 3)
fix5 <- geewa(
  formula = numobese ~ gender,
  id = id,
  data = muscatine,
  family = binomial("logit"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref5$geese$gamma,
  alpha_vector = R5[lower.tri(R5)]
)

check_geepack_fixed_example(
  label = "Ex5 muscatine binomial/logit independence fixed from geepack",
  ref = ref5,
  fixed = fix5,
  max_size = 3
)

# =============================================================================
# Example 6: muscatine - binomial, independence (full model)
# numobese ~ gender + cage + cage2 + gender:cage + gender:cage2
# =============================================================================

ref6 <- geeglm(
  numobese ~ gender + cage + cage2 + gender:cage + gender:cage2,
  id = id,
  waves = occasion,
  data = muscatine,
  family = binomial(),
  corstr = "independence"
)
R6 <- extract_corr_matrix(ref6, 3)
fix6 <- geewa(
  formula = numobese ~ gender + cage + cage2 + gender:cage + gender:cage2,
  id = id,
  data = muscatine,
  family = binomial("logit"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref6$geese$gamma,
  alpha_vector = R6[lower.tri(R6)]
)

check_geepack_fixed_example(
  label = "Ex6 muscatine binomial/logit independence fixed from geepack",
  ref = ref6,
  fixed = fix6,
  max_size = 3
)

# =============================================================================
# Example 7: respiratory - binomial, independence
# =============================================================================

data(respiratory, package = "geepack")
respiratory$center <- factor(respiratory$center)
respiratory$id <- rep(1:111, each = 4)
ref7 <- geeglm(
  outcome ~ center + treat + age + baseline,
  id = id,
  data = respiratory,
  family = binomial(),
  corstr = "independence"
)
R7 <- extract_corr_matrix(ref7, 4)
fix7 <- geewa(
  formula = outcome ~ center + treat + age + baseline,
  id = id,
  data = respiratory,
  family = binomial("logit"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref7$geese$gamma,
  alpha_vector = R7[lower.tri(R7)]
)

check_geepack_fixed_example(
  label = "Ex7 respiratory binomial/logit independence fixed from geepack",
  ref = ref7,
  fixed = fix7,
  max_size = 4
)

# =============================================================================
# Example 8: respiratory - binomial, exchangeable
# =============================================================================

ref8 <- geeglm(
  outcome ~ center + treat + age + baseline,
  id = id,
  data = respiratory,
  family = binomial(),
  corstr = "exchangeable"
)
R8 <- extract_corr_matrix(ref8, 4)
fix8 <- geewa(
  formula = outcome ~ center + treat + age + baseline,
  id = id,
  data = respiratory,
  family = binomial("logit"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref8$geese$gamma,
  alpha_vector = R8[lower.tri(R8)]
)

check_geepack_fixed_example(
  label = "Ex8 respiratory binomial/logit exchangeable fixed from geepack",
  ref = ref8,
  fixed = fix8,
  max_size = 4
)

# =============================================================================
# Example 9: respiratory - binomial, unstructured
# =============================================================================

ref9 <- geeglm(
  outcome ~ center + treat + age + baseline,
  id = id,
  data = respiratory,
  family = binomial(),
  corstr = "unstructured"
)
R9 <- extract_corr_matrix(ref9, 4)
fix9 <- geewa(
  formula = outcome ~ center + treat + age + baseline,
  id = id,
  data = respiratory,
  family = binomial("logit"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref9$geese$gamma,
  alpha_vector = R9[lower.tri(R9)]
)

check_geepack_fixed_example(
  label = "Ex9 respiratory binomial/logit unstructured fixed from geepack",
  ref = ref9,
  fixed = fix9,
  max_size = 4
)

# =============================================================================
# Example 10: respiratory - binomial, ar1
# =============================================================================

ref10 <- geeglm(
  outcome ~ center + treat + age + baseline,
  id = id,
  data = respiratory,
  family = binomial(),
  corstr = "ar1"
)
R10 <- extract_corr_matrix(ref10, 4)
fix10 <- geewa(
  formula = outcome ~ center + treat + age + baseline,
  id = id,
  data = respiratory,
  family = binomial("logit"),
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref10$geese$gamma,
  alpha_vector = R10[lower.tri(R10)]
)

check_geepack_fixed_example(
  label = "Ex10 respiratory binomial/logit ar1 fixed from geepack",
  ref = ref10,
  fixed = fix10,
  max_size = 4
)
