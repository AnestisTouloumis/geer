library(glmtoolbox)
library(testthat)
devtools::load_all()

strip <- function(x) {
  dimnames(x) <- NULL
  x
}

check_glmtoolbox_example <- function(label, ref, fixed, tol = 1e-5) {
  naive_ref <- strip(ref$naive * ref$phi)

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

check_glmtoolbox_criteria <- function(label, ref, fixed, tol = 1e-5) {
  describe(paste(label, "criteria"), {
    it("QIC agrees", {
      expect_equal(
        unname(glmtoolbox::QIC(ref)),
        round(unname(geecriteria(fixed, cov_type = "robust", digits = 3)[1, "QIC"]), 3),
        tolerance = tol
      )
    })
    it("CIC agrees", {
      expect_equal(
        unname(glmtoolbox::CIC(ref, digits = 3)),
        round(unname(geecriteria(fixed, cov_type = "robust", digits = 3)[1, "CIC"]), 3),
        tolerance = tol
      )
    })
    it("RJC agrees", {
      expect_equal(
        unname(glmtoolbox::RJC(ref, digits = 3)),
        round(unname(geecriteria(fixed, cov_type = "robust", digits = 3)[1, "RJC"]), 3),
        tolerance = tol
      )
    })
  })
}

check_confint_example <- function(label, ref, fixed, level = 0.95, tol = 1e-5, ...) {
  ci_ref <- confint(ref, level = level, ...)
  ci_fixed <- confint(fixed, level = level, cov_type = "robust", ...)

  describe(paste(label, "confidence intervals"), {
    it("confidence intervals agree", {
      expect_equal(
        unname(ci_ref),
        unname(ci_fixed),
        tolerance = tol
      )
    })
  })
}

check_glmtoolbox_residuals <- function(label, ref, fixed, tol = 1e-5, ...) {
  res_ref <- residuals(ref, type = "pearson", ...)
  res_fixed <- residuals.geer(fixed, type = "pearson", ...)/sqrt(ref$phi)

  describe(paste(label, "pearson residuals"), {
    it("residuals agree", {
      expect_equal(
        as.vector(unname(res_ref)),
        as.vector(unname(res_fixed)),
        tolerance = tol
      )
    })
  })
}

check_glmtoolbox_fitted <- function(label, ref, fixed, tol = 1e-5, ...) {
  fit_ref <- fitted(ref, ...)
  fit_fixed <- fitted(fixed, ...)

  describe(paste(label, "fitted values"), {
    it("fitted values agree", {
      expect_equal(
        as.vector(unname(fit_ref)),
        as.vector(unname(fit_fixed)),
        tolerance = tol
      )
    })
  })
}

# =============================================================================
# Example 1: Effect of ozone-enriched atmosphere on growth of sitka spruces
# Gamma(log) | AR-M-dependent(1) -> ar1
# =============================================================================

data(spruces)

mod1 <- size ~ poly(days, 4) + treat

ref1 <- glmgee(
  mod1,
  id = tree,
  family = Gamma(log),
  corstr = "AR-M-dependent(1)",
  data = spruces
)

R1 <- ref1$corr

fix1 <- geewa(
  formula = mod1,
  id = tree,
  family = Gamma("log"),
  data = spruces,
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref1$phi,
  alpha_vector = R1[lower.tri(R1)]
)

check_glmtoolbox_example(
  label = "Ex1 spruces Gamma/log ar1",
  ref = ref1,
  fixed = fix1
)

check_glmtoolbox_criteria(
  label = "Ex1 spruces Gamma/log ar1",
  ref = ref1,
  fixed = fix1
)

check_confint_example(
  label = "Ex1 spruces Gamma/log ar1",
  ref = ref1,
  fixed = fix1
)

check_glmtoolbox_residuals(
  label = "Ex1 spruces Gamma/log ar1",
  ref = ref1,
  fixed = fix1
)

check_glmtoolbox_fitted(
  label = "Ex1 spruces Gamma/log ar1",
  ref = ref1,
  fixed = fix1
)


# =============================================================================
# Example 2: Treatment for severe postnatal depression
# binomial(logit) | AR-M-dependent(1) -> ar1
# =============================================================================

data(depression, package = "glmtoolbox")

mod2 <- depressd ~ visit + group

ref2 <- glmgee(
  mod2,
  id = subj,
  family = binomial(logit),
  corstr = "AR-M-dependent(1)",
  data = depression
)

R2 <- ref2$corr

fix2 <- geewa(
  formula = mod2,
  id = subj,
  family = binomial("logit"),
  data = depression,
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref2$phi,
  alpha_vector = R2[lower.tri(R2)]
)

check_glmtoolbox_example(
  label = "Ex2 depression binomial/logit ar1",
  ref = ref2,
  fixed = fix2
)

check_glmtoolbox_criteria(
  label = "Ex2 depression binomial/logit ar1",
  ref = ref2,
  fixed = fix2
)

check_confint_example(
  label = "Ex2 depression binomial/logit ar1",
  ref = ref2,
  fixed = fix2
)

check_glmtoolbox_residuals(
  label = "Ex2 depression binomial/logit ar1",
  ref = ref2,
  fixed = fix2
)

check_glmtoolbox_fitted(
  label = "Ex2 depression binomial/logit ar1",
  ref = ref2,
  fixed = fix2
)

# =============================================================================
# Example 3: Treatment for severe postnatal depression (2)
# gaussian | AR-M-dependent(1) -> ar1
# =============================================================================

mod3 <- dep ~ visit * group

ref3 <- glmgee(
  mod3,
  id = subj,
  family = gaussian,
  corstr = "AR-M-dependent(1)",
  data = depression
)

R3 <- ref3$corr

fix3 <- geewa(
  formula = mod3,
  id = subj,
  family = gaussian,
  data = depression,
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref3$phi,
  alpha_vector = R3[lower.tri(R3)]
)

check_glmtoolbox_example(
  label = "Ex3 depression gaussian ar1",
  ref = ref3,
  fixed = fix3
)

check_glmtoolbox_criteria(
  label = "Ex3 depression gaussian ar1",
  ref = ref3,
  fixed = fix3
)

check_confint_example(
  label = "Ex3 depression gaussian ar1",
  ref = ref3,
  fixed = fix3
)

check_glmtoolbox_residuals(
  label = "Ex3 depression gaussian ar1",
  ref = ref3,
  fixed = fix3
)

check_glmtoolbox_fitted(
  label = "Ex3 depression gaussian ar1",
  ref = ref3,
  fixed = fix3
)

# =============================================================================
# Example 4: Dental Clinical Trial
# binomial(log) | Exchangeable
# =============================================================================

data(rinse, package = "glmtoolbox")
rinse$id_new <- match(rinse$subject, unique(rinse$subject))

mod4 <- score / 3.6 ~ rinse * time


ref4 <- glmgee(
  mod4,
  id = id_new,
  family = binomial(log),
  corstr = "Exchangeable",
  data = rinse
)

R4 <- ref4$corr

fix4 <- geewa(
  formula = mod4,
  id = id_new,
  family = binomial("log"),
  data = rinse,
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref4$phi,
  alpha_vector = R4[lower.tri(R4)]
)

check_glmtoolbox_example(
  label = "Ex4 rinse binomial/log exchangeable",
  ref = ref4,
  fixed = fix4
)

check_glmtoolbox_criteria(
  label = "Ex4 rinse binomial/log exchangeable",
  ref = ref4,
  fixed = fix4
)

check_confint_example(
  label = "Ex4 rinse binomial/log exchangeable",
  ref = ref4,
  fixed = fix4,
  tol = 1e-4
)

check_glmtoolbox_residuals(
  label = "Ex4 rinse binomial/log exchangeable",
  ref = ref4,
  fixed = fix4
)

check_glmtoolbox_fitted(
  label = "Ex4 rinse binomial/log exchangeable",
  ref = ref4,
  fixed = fix4
)

# =============================================================================
# Example 5: Shoulder Pain after Laparoscopic Cholecystectomy
# binomial(logit) | Stationary-M-dependent(2)
# =============================================================================

data(cholecystectomy, package = "glmtoolbox")

mod5 <- pain2 ~ treatment + age + time

ref5 <- glmgee(
  mod5,
  id = id,
  family = binomial(logit),
  corstr = "Stationary-M-dependent(2)",
  data = cholecystectomy
)

R5 <- ref5$corr

fix5 <- geewa(
  formula = mod5,
  id = id,
  family = binomial("logit"),
  data = cholecystectomy,
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref5$phi,
  alpha_vector = R5[lower.tri(R5)]
)

check_glmtoolbox_example(
  label = "Ex5 cholecystectomy binomial/logit m-dep(2)",
  ref = ref5,
  fixed = fix5
)

check_glmtoolbox_criteria(
  label = "Ex5 cholecystectomy binomial/logit m-dep(2)",
  ref = ref5,
  fixed = fix5
)

check_confint_example(
  label = "Ex5 cholecystectomy binomial/logit m-dep(2)",
  ref = ref5,
  fixed = fix5
)

check_glmtoolbox_residuals(
  label = "Ex5 cholecystectomy binomial/logit m-dep(2)",
  ref = ref5,
  fixed = fix5
)

check_glmtoolbox_fitted(
  label = "Ex5 cholecystectomy binomial/logit m-dep(2)",
  ref = ref5,
  fixed = fix5
)

# =============================================================================
# Example 6: Guidelines for Urinary Incontinence Discussion and Evaluation
# binomial(logit) | Exchangeable
# =============================================================================

data(GUIDE, package = "glmtoolbox")
GUIDE$id_new <- match(GUIDE$practice, unique(GUIDE$practice))

mod6 <- bothered ~ gender + age + dayacc + severe + toilet

ref6 <- glmgee(
  mod6,
  id = id_new,
  family = binomial(logit),
  corstr = "Exchangeable",
  data = GUIDE
)

R6 <- ref6$corr

fix6 <- geewa(
  formula = mod6,
  id = id_new,
  family = binomial("logit"),
  data = GUIDE,
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref6$phi,
  alpha_vector = R6[lower.tri(R6)]
)

check_glmtoolbox_example(
  label = "Ex6 GUIDE binomial/logit exchangeable",
  ref = ref6,
  fixed = fix6
)

check_glmtoolbox_criteria(
  label = "Ex6 GUIDE binomial/logit exchangeable",
  ref = ref6,
  fixed = fix6
)

check_confint_example(
  label = "Ex6 GUIDE binomial/logit exchangeable",
  ref = ref6,
  fixed = fix6
)

check_glmtoolbox_residuals(
  label = "Ex6 GUIDE binomial/logit exchangeable",
  ref = ref6,
  fixed = fix6
)

check_glmtoolbox_fitted(
  label = "Ex6 GUIDE binomial/logit exchangeable",
  ref = ref6,
  fixed = fix6
)

# =============================================================================
# Example 7: Tests of Auditory Perception in Children with OME
# binomial(cloglog) | Exchangeable
# =============================================================================

OME <- MASS::OME
OME$id_new <- match(OME$ID, unique(OME$ID))
mod7 <- cbind(Correct, Trials - Correct) ~ Loud + Age + OME

ref7 <- glmgee(
  mod7,
  id = id_new,
  family = binomial(cloglog),
  corstr = "Exchangeable",
  data = OME
)

R7 <- ref7$corr

fix7 <- geewa(
  formula = mod7,
  id = id_new,
  family = binomial("cloglog"),
  data = OME,
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref7$phi,
  alpha_vector = R7[lower.tri(R7)]
)

check_glmtoolbox_example(
  label = "Ex7 OME binomial/cloglog exchangeable",
  ref = ref7,
  fixed = fix7
)

check_glmtoolbox_criteria(
  label = "Ex7 OME binomial/cloglog exchangeable",
  ref = ref7,
  fixed = fix7
)

check_confint_example(
  label = "Ex7 OME binomial/cloglog exchangeable",
  ref = ref7,
  fixed = fix7
)

check_glmtoolbox_residuals(
  label = "Ex7 OME binomial/cloglog exchangeable",
  ref = ref7,
  fixed = fix7
)

check_glmtoolbox_fitted(
  label = "Ex7 OME binomial/cloglog exchangeable",
  ref = ref7,
  fixed = fix7
)

# =============================================================================
# Example 8: Epileptic seizures
# poisson(log) | Exchangeable
# =============================================================================

data(Seizures, package = "glmtoolbox")

Seizures2 <- within(Seizures, time4 <- ifelse(time == 4, 1, 0))
mod8 <- seizures ~ log(age) + time4 + log(base / 4) * treatment

ref8 <- glmgee(
  mod8,
  id = id,
  family = poisson(log),
  corstr = "Exchangeable",
  data = Seizures2
)

R8 <- ref8$corr

fix8 <- geewa(
  formula = mod8,
  id = id,
  family = poisson("log"),
  data = Seizures2,
  corstr = "fixed",
  phi_fixed = TRUE,
  phi_value = ref8$phi,
  alpha_vector = R8[lower.tri(R8)]
)

check_glmtoolbox_example(
  label = "Ex8 Seizures poisson/log exchangeable",
  ref = ref8,
  fixed = fix8
)

check_glmtoolbox_criteria(
  label = "Ex8 Seizures poisson/log exchangeable",
  ref = ref8,
  fixed = fix8
)

check_confint_example(
  label = "Ex8 Seizures poisson/log exchangeable",
  ref = ref8,
  fixed = fix8
)

check_glmtoolbox_residuals(
  label = "Ex8 Seizures poisson/log exchangeable",
  ref = ref8,
  fixed = fix8
)

check_glmtoolbox_fitted(
  label = "Ex8 Seizures poisson/log exchangeable",
  ref = ref8,
  fixed = fix8
)
