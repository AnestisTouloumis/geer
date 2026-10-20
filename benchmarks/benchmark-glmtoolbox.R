# =============================================================================
# Benchmark: glmtoolbox (glmgee) versus geer::geewa
# =============================================================================

library(glmtoolbox)
library(testthat)
bench_dir <- if (file.exists("benchmark-helpers.R")) "." else "benchmarks"
devtools::load_all(file.path(bench_dir, ".."))
source(file.path(bench_dir, "benchmark-helpers.R"))

glmgee_fit <- function(formula, id, family, data, corstr = "Independence") {
  eval(as.call(list(quote(glmgee), formula = formula, family = family,
                    id = as.name(id), data = data, corstr = corstr,
                    toler = 1e-10, maxit = 200)))
}

geer_cluster_ids <- function(fit, id) {
  used <- id[fit$frame_rows]
  unique(as.character(used[fit$row_order]))
}

fixed_from_glmgee <- function(args, ref) {
  geewa_fixed(args, ref$corr, ref$phi)
}

glmgee_m <- function(ref) {
  corstr <- ref$corstr
  if (identical(as.character(corstr), "Independence")) return(0)
  if (identical(as.character(corstr), "Exchangeable")) return(1)
  attr(corstr, "M")
}

as_ref_glmgee <- function(ref, fit, id) {
  mahal <- num(residuals(ref, type = "mahalanobis"))
  data_order <- unique(as.character(id[fit$frame_rows]))
  mahal <- mahal[match(geer_cluster_ids(fit, id), data_order)]
  list(
    coefficients = coef(ref),
    naive = vcov(ref, type = "model"),
    robust = vcov(ref, type = "robust"),
    df_adjusted = vcov(ref, type = "df-adjusted"),
    fitted = ref$fitted.values,
    linear_predictors = ref$linear.predictors,
    residuals_pearson = residuals(ref, type = "pearson"),
    residuals_deviance = residuals(ref, type = "deviance"),
    residuals_mahalanobis = mahal,
    confint = list(cov_type = "robust", level = 0.95,
                   ci = confint(ref, verbose = FALSE, digits = 15))
  )
}

check_criteria_glmgee <- function(label, ref, fit) {
  m <- glmgee_m(ref)
  K <- fit$clusters_no
  g <- geer_criteria(fit, c("AGPC", "SGPC"))
  unbalanced <- "geer requires a balanced design for GHYC and PAC"
  balanced <- all(table(fit$id) == length(unique(fit$repeated)))
  check_criteria(
    paste(label, "criteria"),
    c(QIC = -2 * ref$logLik + 2 * ref$CIC,
      QICu = -2 * ref$logLik + 2 * length(coef(ref)),
      CIC = CIC(ref, verbose = FALSE, digits = 15),
      RJC = RJC(ref, verbose = FALSE, digits = 15),
      GHYC = GHYC(ref, verbose = FALSE, digits = 15),
      PAC = PAC(ref, verbose = FALSE, digits = 15)),
    fit,
    known = if (!balanced) list(GHYC = unbalanced, PAC = unbalanced)
  )
  describe(paste(label, "AGPC and SGPC"), {
    it("AGPC(glmgee) = AGPC(geer) + 2 m", {
      expect_agree(AGPC(ref, verbose = FALSE, digits = 15),
                   g[["AGPC"]] + 2 * m)
    })
    it("SGPC(glmgee) = SGPC(geer) + log(K) m", {
      expect_agree(SGPC(ref, verbose = FALSE, digits = 15),
                   g[["SGPC"]] + log(K) * m)
    })
  })
}

check_glmgee <- function(label, ref, fit, id) {
  check_fit(label, as_ref_glmgee(ref, fit, id), fit)
  check_criteria_glmgee(label, ref, fit)
}

check_predict_glmgee <- function(label, ref, fit, newdata) {
  p_ref <- predict(ref, newdata = newdata, type = "response", se.fit = TRUE)
  p_fit <- predict(fit, newdata = newdata, type = "response", se.fit = TRUE,
                   cov_type = "robust")
  describe(paste(label, "predict()"), {
    it("predictions agree", {
      expect_agree(p_ref[, "fit"], p_fit$fit)
    })
    it("standard errors agree", {
      expect_agree(p_ref[, "se.fit"], p_fit$se.fit)
    })
  })
}

anova_relation <- function(refs, formulas, args, test) {
  vapply(seq_along(refs)[-1L], function(i) {
    ref_i <- if (test == "wald") refs[[i]] else refs[[i - 1L]]
    small <- fixed_from_glmgee(c(list(formula = formulas[[i - 1L]]), args),
                               ref_i)
    large <- fixed_from_glmgee(c(list(formula = formulas[[i]]), args), ref_i)
    anova(small, large, test = test, cov_type = "robust")$Chi[2L]
  }, numeric(1))
}

sequential_formulas <- function(formula) {
  response <- deparse(formula[[2L]])
  labels <- attr(terms(formula), "term.labels")
  c(list(reformulate("1", response = response)),
    lapply(seq_along(labels), function(i) {
      reformulate(labels[seq_len(i)], response = response)
    }))
}


# =============================================================================
# Example 1: spruces, Gamma(log)
# =============================================================================

data("spruces", package = "glmtoolbox")
mod1 <- size ~ poly(days, 4) + treat
args_sp <- list(id = quote(tree), data = spruces,
                family = Gamma(link = "log"))

ref_sp <- list(
  independence = glmgee_fit(mod1, "tree", Gamma(log), spruces),
  "AR-M-dependent(1)" = glmgee_fit(mod1, "tree", Gamma(log), spruces,
                                   "AR-M-dependent(1)"),
  exchangeable = glmgee_fit(mod1, "tree", Gamma(log), spruces,
                            "Exchangeable")
)
fit_sp <- lapply(ref_sp, function(ref) {
  fixed_from_glmgee(c(list(formula = mod1), args_sp), ref)
})
for (cs in names(ref_sp)) {
  check_glmgee(paste("Ex1 spruces Gamma/log", cs), ref_sp[[cs]],
               fit_sp[[cs]], spruces$tree)
}
not_comparable(
  "Ex1 spruces Stationary-M-dependent(2) (criteria examples, fit3)",
  paste("the estimated 2-dependent correlation matrix is not positive",
        "definite, so it cannot be used as a fixed working correlation")
)

newdata_sp <- data.frame(days = c(556, 556),
                         treat = as.factor(c("normal", "ozone-enriched")))
check_predict_glmgee("Ex1 spruces AR-M-dependent(1)",
                     ref_sp[["AR-M-dependent(1)"]],
                     fit_sp[["AR-M-dependent(1)"]], newdata_sp)

formulas_a <- list(size ~ poly(days, 4),
                   size ~ poly(days, 4) + treat,
                   size ~ poly(days, 4) + treat + poly(days, 4):treat)
refs_a <- lapply(formulas_a, glmgee_fit, id = "tree", family = Gamma(log),
                 data = spruces, corstr = "AR-M-dependent")
seq_a <- sequential_formulas(formulas_a[[3L]])
refs_seq_a <- lapply(seq_a, glmgee_fit, id = "tree", family = Gamma(log),
                     data = spruces, corstr = "AR-M-dependent")

describe("Ex1 spruces anova() Wald tests", {
  it("nested comparisons agree", {
    a_ref <- anova(refs_a[[1L]], refs_a[[2L]], refs_a[[3L]], test = "wald",
                   verbose = FALSE)
    expect_agree(a_ref[, 1L],
                 anova_relation(refs_a, formulas_a, args_sp, "wald"))
  })
  it("sequential tests agree", {
    a_ref <- anova(refs_a[[3L]], test = "wald", verbose = FALSE)
    expect_agree(a_ref[, 1L],
                 anova_relation(refs_seq_a, seq_a, args_sp, "wald"))
  })
})


# =============================================================================
# Example 2: depression, binomial(logit)
# =============================================================================

data("depression", package = "glmtoolbox")
mod2 <- depressd ~ visit + group
args_d2 <- list(id = quote(subj), data = depression,
                family = binomial(link = "logit"))

ref_d2 <- list(
  independence = glmgee_fit(mod2, "subj", binomial(logit), depression),
  "AR-M-dependent(1)" = glmgee_fit(mod2, "subj", binomial(logit), depression,
                                   "AR-M-dependent(1)"),
  "Stationary-M-dependent(2)" = glmgee_fit(mod2, "subj", binomial(logit),
                                           depression,
                                           "Stationary-M-dependent(2)"),
  exchangeable = glmgee_fit(mod2, "subj", binomial(logit), depression,
                            "Exchangeable")
)
fit_d2 <- lapply(ref_d2, function(ref) {
  fixed_from_glmgee(c(list(formula = mod2), args_d2), ref)
})
for (cs in names(ref_d2)) {
  check_glmgee(paste("Ex2 depression binomial/logit", cs), ref_d2[[cs]],
               fit_d2[[cs]], depression$subj)
}

newdata_d <- data.frame(visit = c(6, 6),
                        group = as.factor(c("placebo", "oestrogen")))
check_predict_glmgee("Ex2 depression AR-M-dependent(1)",
                     ref_d2[["AR-M-dependent(1)"]],
                     fit_d2[["AR-M-dependent(1)"]], newdata_d)

formulas_b <- list(depressd ~ group,
                   depressd ~ group + visit,
                   depressd ~ group + visit + group:visit)
refs_b <- lapply(formulas_b, glmgee_fit, id = "subj",
                 family = binomial(logit), data = depression,
                 corstr = "AR-M-dependent")
seq_b <- sequential_formulas(formulas_b[[3L]])
refs_seq_b <- lapply(seq_b, glmgee_fit, id = "subj",
                     family = binomial(logit), data = depression,
                     corstr = "AR-M-dependent")

describe("Ex2 depression anova() score tests", {
  it("nested comparisons agree", {
    a_ref <- anova(refs_b[[1L]], refs_b[[2L]], refs_b[[3L]], test = "score",
                   verbose = FALSE)
    expect_agree(a_ref[, 1L],
                 anova_relation(refs_b, formulas_b, args_d2, "score"))
  })
  it("sequential tests agree", {
    a_ref <- anova(refs_b[[3L]], test = "score", verbose = FALSE)
    expect_agree(a_ref[, 1L],
                 anova_relation(refs_seq_b, seq_b, args_d2, "score"))
  })
})


# =============================================================================
# Example 3: depression, gaussian(identity)
# =============================================================================

mod3 <- dep ~ visit * group
args_d3 <- list(id = quote(subj), data = depression,
                family = gaussian(link = "identity"))

ref_d3 <- list(
  independence = glmgee_fit(mod3, "subj", gaussian(identity), depression),
  "AR-M-dependent(1)" = glmgee_fit(mod3, "subj", gaussian(identity),
                                   depression, "AR-M-dependent(1)"),
  exchangeable = glmgee_fit(mod3, "subj", gaussian(identity), depression,
                            "Exchangeable")
)
fit_d3 <- lapply(ref_d3, function(ref) {
  fixed_from_glmgee(c(list(formula = mod3), args_d3), ref)
})
for (cs in names(ref_d3)) {
  check_glmgee(paste("Ex3 depression gaussian", cs), ref_d3[[cs]],
               fit_d3[[cs]], depression$subj)
}
check_predict_glmgee("Ex3 depression gaussian AR-M-dependent(1)",
                     ref_d3[["AR-M-dependent(1)"]],
                     fit_d3[["AR-M-dependent(1)"]], newdata_d)


# =============================================================================
# Example 4: rinse, binomial(log), exchangeable
# =============================================================================

data("rinse", package = "glmtoolbox")
ref4 <- glmgee_fit(score / 3.6 ~ rinse * time, "subject", binomial(log),
                   rinse, "Exchangeable")
fit4 <- fixed_from_glmgee(list(formula = score / 3.6 ~ rinse * time,
                               id = quote(subject), data = rinse,
                               family = binomial(link = "log")), ref4)
check_glmgee("Ex4 rinse binomial/log exchangeable", ref4, fit4,
             rinse$subject)


# =============================================================================
# Example 5: cholecystectomy, binomial(logit), Stationary-M-dependent(2)
# =============================================================================

data("cholecystectomy", package = "glmtoolbox")
ref5 <- glmgee_fit(pain2 ~ treatment + age + time, "id", binomial(logit),
                   cholecystectomy, "Stationary-M-dependent(2)")
fit5 <- fixed_from_glmgee(list(formula = pain2 ~ treatment + age + time,
                               id = quote(id), data = cholecystectomy,
                               family = binomial(link = "logit")), ref5)
check_glmgee("Ex5 cholecystectomy binomial/logit m-dependent(2)", ref5,
             fit5, cholecystectomy$id)


# =============================================================================
# Example 6: GUIDE, binomial(logit), exchangeable (also ?GUIDE, ?leverage)
# =============================================================================

data("GUIDE", package = "glmtoolbox")
mod6 <- bothered ~ gender + age + dayacc + severe + toilet
ref6 <- glmgee_fit(mod6, "practice", binomial(logit), GUIDE, "Exchangeable")
fit6 <- fixed_from_glmgee(list(formula = mod6, id = quote(practice),
                               data = GUIDE,
                               family = binomial(link = "logit")), ref6)
check_glmgee("Ex6 GUIDE binomial/logit exchangeable", ref6, fit6,
             GUIDE$practice)


# =============================================================================
# Example 7: OME (MASS), binomial(cloglog), exchangeable (also ?leverage)
# =============================================================================

OME <- MASS::OME
mod7 <- cbind(Correct, Trials - Correct) ~ Loud + Age + OME
ref7 <- glmgee_fit(mod7, "ID", binomial(cloglog), OME, "Exchangeable")
fit7 <- fixed_from_glmgee(list(formula = mod7, id = quote(ID), data = OME,
                               family = binomial(link = "cloglog")), ref7)
check_glmgee("Ex7 OME binomial/cloglog exchangeable", ref7, fit7, OME$ID)


# =============================================================================
# Example 8: Seizures, poisson(log), exchangeable
# =============================================================================

data("Seizures", package = "glmtoolbox")
Seizures2 <- within(Seizures, time4 <- ifelse(time == 4, 1, 0))
mod8 <- seizures ~ log(age) + time4 + log(base / 4) * treatment
ref8 <- glmgee_fit(mod8, "id", poisson(log), Seizures2, "Exchangeable")
fit8 <- fixed_from_glmgee(list(formula = mod8, id = quote(id),
                               data = Seizures2,
                               family = poisson(link = "log")), ref8)
check_glmgee("Ex8 Seizures poisson/log exchangeable", ref8, fit8,
             Seizures2$id)


# =============================================================================
# stepCriterion() versus step_criterion()
# =============================================================================

sorted_term <- function(term) {
  vapply(
    strsplit(term, ":", fixed = TRUE),
    function(x) paste(sort(x), collapse = ":"),
    character(1)
  )
}

glmgee_steps <- function(trace_output) {
  moves <- grep("^Step [0-9]+ : [+-] ", trace_output, value = TRUE)
  moves <- sub("^Step [0-9]+ : ", "", trimws(moves))
  paste(substr(moves, 1L, 2L), sorted_term(substring(moves, 3L)))
}

geer_steps <- function(fit) {
  moves <- rownames(fit$anova)[-1L]
  moves <- sub("^Step [0-9]+: ", "", moves)
  paste(substr(moves, 1L, 2L), sorted_term(substring(moves, 3L)))
}

glmgee_criterion <- function(ref, criterion) {
  as.numeric(switch(
    criterion,
    QIC = -2 * ref$logLik + 2 * ref$CIC,
    QICu = -2 * ref$logLik + 2 * length(coef(ref)),
    AGPC = AGPC(ref, verbose = FALSE, digits = 15),
    SGPC = SGPC(ref, verbose = FALSE, digits = 15)
  ))
}

check_step_criterion <- function(label, start_formula, lower, upper, id,
                                 family, data, direction) {
  for (criterion in c("QIC", "QICu", "AGPC", "SGPC")) {
    ref_start <- glmgee_fit(start_formula, id, family, data)
    trace_output <- capture.output(
      ref <- stepCriterion(
        ref_start,
        criterion = tolower(criterion),
        direction = direction,
        scope = list(lower = lower, upper = upper),
        force.in = lower,
        trace = TRUE
      )
    )
    ref_final <- glmgee_fit(
      as.formula(paste(deparse(start_formula[[2L]]), ref$final)),
      id, family, data
    )
    fit_start <- geewa_native(list(
      formula = start_formula, id = as.name(id), data = data,
      family = family, corstr = "independence"
    ))
    fit <- step_criterion(
      fit_start,
      scope = list(lower = lower, upper = upper),
      direction = "both",
      criterion = criterion
    )
    describe(paste(label, criterion), {
      it("the sequence of moves agrees", {
        expect_identical(geer_steps(fit), glmgee_steps(trace_output))
      })
      it("the selected model agrees", {
        expect_setequal(
          sorted_term(attr(terms(formula(fit)), "term.labels")),
          sorted_term(attr(terms(formula(ref_final)), "term.labels"))
        )
      })
      it("the criterion of the initial model agrees", {
        expect_agree(glmgee_criterion(ref_start, criterion),
                     fit$anova[[criterion]][1L])
      })
      it("the criterion of the selected model agrees", {
        expect_agree(glmgee_criterion(ref_final, criterion),
                     fit$anova[[criterion]][nrow(fit$anova)])
      })
      it("the criterion of every candidate agrees", {
        for (candidate_table in fit$candidates) {
          current <- as.formula(attr(candidate_table, "heading")[3L])
          changes <- setdiff(rownames(candidate_table), "<none>")
          for (change in changes) {
            ref_candidate <- glmgee_fit(
              update(current, as.formula(paste(". ~ .", change))),
              id, family, data
            )
            expect_agree(glmgee_criterion(ref_candidate, criterion),
                         candidate_table[change, criterion])
          }
        }
      })
    })
  }
}

data("depression", package = "glmtoolbox")
check_step_criterion(
  "?stepCriterion Ex2 depression binomial/probit forward",
  start_formula = depressd ~ group,
  lower = ~group,
  upper = depressd ~ visit * group,
  id = "subj",
  family = binomial(probit),
  data = depression,
  direction = "forward"
)
check_step_criterion(
  "?stepCriterion Ex3 depression gaussian/identity forward",
  start_formula = dep ~ group,
  lower = ~group,
  upper = dep ~ visit * group,
  id = "subj",
  family = gaussian(identity),
  data = depression,
  direction = "forward"
)
data("spruces", package = "glmtoolbox")
check_step_criterion(
  "?stepCriterion Ex1 spruces Gamma/log forward",
  start_formula = size ~ treat,
  lower = ~treat,
  upper = size ~ poly(days, 4) * treat,
  id = "tree",
  family = Gamma(log),
  data = spruces,
  direction = "forward"
)
not_comparable("stepCriterion(direction = \"backward\") with force.in",
               paste("glmtoolbox stops with an error when no candidate is",
                     "left at the final model"))


# =============================================================================
# Examples without a geer counterpart
# =============================================================================

not_comparable("stepCriterion(criterion = \"p-value\")",
               paste("glmtoolbox mixes additions and deletions with its own",
                     "hybrid rule and Wilson-Hilferty p-values, which step_p()",
                     "does not reproduce"))
not_comparable("cooks.distance(), dfbeta(), leverage(), localInfluence()",
               "geer has no influence diagnostics")
not_comparable("estequa()", "geer does not export the estimating equations")
not_comparable("wglmgee() (amenorrhea)",
               "geer has no weighted GEE for dropout")
