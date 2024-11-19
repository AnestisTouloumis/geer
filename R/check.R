library("tidyverse")

devtools::load_all()
obesity <-
  obesity |>
  mutate(
    y = if_else(obesity == "Lean", 1, 0)
  )

fitted_model_geepack <-
  glmtoolbox::glmgee(
    formula = y ~ bacteroides + ancestry + age + zygosity,
    id = fid,
    data = obesity,
    family = binomial(link = "logit"),
    corstr = "independence"
  )

fitted_model_gee <-
  geewa(
    formula = y ~ bacteroides + ancestry + age + zygosity,
    id = fid,
    data = obesity,
    correlation_structure = "fixed",
    family = binomial(link = "logit"),
    alpha_vector = fitted_model_geepack$corr[lower.tri(fitted_model_geepack$corr)],
    method = "gee",
    phi_value = fitted_model_geepack$phi,
    phi_fixed = TRUE
  )
fitted_model_brgee_robust <-
  update(fitted_model_gee, method = "brgee_robust")
fitted_model_bcgee_robust <-
  update(fitted_model_gee, method = "bcgee_robust")
fitted_model_brgee_naive <-
  update(fitted_model_gee, method = "brgee_naive")
fitted_model_bcgee_naive <-
  update(fitted_model_gee, method = "bcgee_naive")
gee_criteria(fitted_model_gee,
             fitted_model_brgee_robust,
             fitted_model_bcgee_robust,
             fitted_model_brgee_naive,
             fitted_model_bcgee_naive,
             cov_type = "robust", digits = 3
)
glmtoolbox::QIC(fitted_model_geepack)
glmtoolbox::CIC(fitted_model_geepack)
glmtoolbox::RJC(fitted_model_geepack)
gee_criteria(fitted_model_gee, digits  = 5)




gee_criteria(fitted_model_gee,
             fitted_model_brgee_robust,
             fitted_model_bcgee_robust,
             fitted_model_brgee_naive,
             fitted_model_bcgee_naive,
             cov_type = "bias-corrected", digits = 4
)



wheeze <- geessbin::wheeze

fitted_model_geepack <-
  glmtoolbox::glmgee(
    formula = Wheeze ~ City + factor(Age),
    data = wheeze,
    id = ID,
    corstr = "exchangeable",
    family = binomial(link = "logit")
  )


fitted_model_gee <- geewa(formula = Wheeze ~ City + factor(Age),
                          data = wheeze,
                          id = ID,
                          correlation_structure = "fixed",
                          repeated = Age,
                          family = binomial(link = "logit"),
                          method = "gee",
                          alpha_vector =
                            fitted_model_geepack$corr[lower.tri(fitted_model_geepack$corr)],
                          phi_fixed = TRUE,
                          phi_value = fitted_model_geepack$phi
)

testthat::expect_equal(fitted(fitted_model_gee), c(fitted(fitted_model_geepack)))

glmtoolbox::QIC(fitted_model_geepack)
glmtoolbox::CIC(fitted_model_geepack)
glmtoolbox::RJC(fitted_model_geepack)

gee_criteria(fitted_model_gee, digits  = 5, cov_type = "robust")
(glmtoolbox::AGPC(fitted_model_geepack) - 2 * 6 - fitted_model_gee$obs_no * log(2 * pi))/2
