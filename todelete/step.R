data("cerebrovascular")
devtools::load_all()
object0 <- geewa(formula = ecg ~ 1,
                 id = id,
                 data = cerebrovascular,
                 family = binomial(link = "logit"),
                 corstr = "independence",
                 method = "pgee_jeffreys",
                 phi_fixed = TRUE)
object1 <- update(object0, formula = ecg ~ period + treatment)
object2 <- update(object0, formula = ecg ~ period*treatment)





step_p(object2,
                       test = "wald",
                       cov_type = "bias-corrected",
                       p_remove = 0.05,
                       direction = "backward")$anova

x1 <- step_p(object0,
                       scope = list(upper = ~ period + treatment + period:treatment),
                       test = "wald",
                       cov_type = "bias-corrected",
                       p_enter = 0.05,
                       direction = "forward")$anova

step_p(object0,
                       scope = list(lower= ~ 1,
                                    upper = ~ period + treatment + period:treatment),
                       test = "wald",
                       cov_type = "bias-corrected",
                       p_enter = 0.05,
                       p_remove = 0.05,
                       direction = "both")$anova

step_p(object2,
                       scope = list(lower= ~ 1,
                                    upper = ~ period + treatment + period:treatment),
                       test = "wald",
                       cov_type = "bias-corrected",
                       p_enter = 0.05,
                       p_remove = 0.05,
                       direction = "both")$anova



step_p(object2,
                       scope = list(upper = ~ period + treatment + period:treatment),
                       test = "wald",
                       p_enter = 0.05,
                       p_remove = 0.05,
                       direction = "both")$anova




object2_glm <- glm(formula = ecg ~ period*treatment,
                   data = cerebrovascular,
                   family = binomial(link = "logit"))
step(object2_glm, direction = "backward", trace = 0)$anova


step(object2_glm,
     scope = list(lower = ~ 1), direction = "both", trace = 0)$anova
