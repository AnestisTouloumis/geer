library("tidyverse")

devtools::load_all()
obesity <-
  obesity |>
  mutate(
    y = if_else(obesity == "Lean", 1, 0)
  )


fitted_model_gee <-
  geewa_binary(
    formula = y ~ bacteroides + ancestry + age + zygosity,
    id = fid,
    data = obesity,
    link = "logit",
    or_structure = "independence",
    method = "gee"
  )
fitted_model_brgee_robust <-
  update(fitted_model_gee, method = "brgee_robust")
fitted_model_bcgee_robust <-
  update(fitted_model_gee, method = "bcgee_robust")
fitted_model_brgee_naive <-
  update(fitted_model_gee, method = "brgee_naive")
fitted_model_bcgee_naive <-
  update(fitted_model_gee, method = "bcgee_naive")
gee_criteria(fitted_model_brgee_robust, cov_type = "robust")
gee_criteria(fitted_model_brgee_robust,
             fitted_model_gee,
             cov_type = "bias-corrected")
gee_criteria(fitted_model_gee,
             fitted_model_brgee_robust,
             fitted_model_bcgee_robust,
             fitted_model_brgee_naive,
             fitted_model_bcgee_naive
             )

fitted_model_geepack <-
  glmtoolbox::glmgee(
    formula = y ~ bacteroides + ancestry + age + zygosity,
    id = fid,
    data = obesity,
    family = binomial(link = "logit"),
    corstr = "Independence"
    )
fitted_model_gee <-
  geewa(
    formula = y ~ bacteroides + ancestry + age + zygosity,
    id = fid,
    data = obesity,
    family = binomial(link = "logit"),
    correlation_structure = "independence",
    method = "gee",
    phi_value = fitted_model_geepack$phi,
    phi_fixed = TRUE
  )
glmtoolbox::QIC(fitted_model_geepack)
glmtoolbox::CIC(fitted_model_geepack)
glmtoolbox::RJC(fitted_model_geepack)
gee_criteria(fitted_model_gee)

glmtoolbox::GHYC(fitted_model_geepack)
glmtoolbox::AGPC(fitted_model_geepack)
glmtoolbox::SGPC(fitted_model_geepack)
fitted_model_geepack$logLik

fitted_glm <-
  glm(
    formula = y ~ bacteroides + ancestry + age + zygosity,
    data = obesity,
    family = binomial(link = "logit")
  )
