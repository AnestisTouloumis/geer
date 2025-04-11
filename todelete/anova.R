data("cerebrovascular")
devtools::load_all()
object0 <- geewa(formula = ecg ~ 1,
                              id = id,
                              data = cerebrovascular,
                              family = binomial(link = "logit"),
                              corstr = "independence",
                              method = "pgee_jeffreys",
                 phi_fixed = TRUE)
object1 <- update(object0, formula = ecg ~ period)
object2 <- update(object0, formula = ecg ~ period*treatment)


cov_types <- c("robust", "naive", "bias-corrected", "df-adjusted")
for(cov_type in cov_types) {
anova(object2, test = "wald", cov_type = cov_type) |> print()
anova(object0, object2, test = "wald", cov_type = cov_type) |> print()
anova(object2, object0, object1, test = "wald", cov_type = cov_type) |> print()
}

cov_types <- c("robust", "naive", "bias-corrected", "df-adjusted")
for(cov_type in cov_types) {
  anova(object1, test = "score", cov_type = cov_type) |> print()
  anova(object0, object1, test = "score", cov_type = cov_type) |> print()
  anova(object0, object1, object2, test = "score", cov_type = cov_type) |> print()
}

pmethods <- c("rao-scott", "satterthwaite")
for(pmethod in pmethods) {
for(cov_type in cov_types) {
  anova(object2, test = "working-wald", cov_type = cov_type, pmethod = pmethod) |> print()
  anova(object0, object2, test = "working-wald", cov_type = cov_type, pmethod = pmethod) |> print()
  anova(object2, object0, object1, test = "working-wald", cov_type = cov_type, pmethod = pmethod) |> print()
  }
}


for(pmethod in pmethods) {
  for(cov_type in cov_types) {
    anova(object2, test = "working-score", cov_type = cov_type, pmethod = pmethod) |> print()
    anova(object0, object2, test = "working-score", cov_type = cov_type, pmethod = pmethod) |> print()
    anova(object2, object0, object1, test = "working-score", cov_type = cov_type, pmethod = pmethod) |> print()
  }
}


for(pmethod in pmethods) {
  for(cov_type in cov_types) {
    anova(object2, test = "working-lrt", cov_type = cov_type, pmethod = pmethod) |> print()
    anova(object0, object2, test = "working-lrt", cov_type = cov_type, pmethod = pmethod) |> print()
    anova(object2, object0, object1, test = "working-lrt", cov_type = cov_type, pmethod = pmethod) |> print()
  }
}


object0 <- geewa_binary(formula = ecg ~ 1,
                 id = id,
                 data = cerebrovascular,
                 link = "logit",
                 orstr = "independence",
                 method = "gee")
object1 <- update(object1, formula = ecg ~ period)
object2 <- update(object1, formula = ecg ~ period * treatment)

cov_types <- c("robust", "naive", "bias-corrected", "df-adjusted")
for(cov_type in cov_types) {
  anova(object2, test = "wald", cov_type = cov_type) |> print()
  anova(object0, object2, test = "wald", cov_type = cov_type) |> print()
  anova(object2, object0, object1, test = "wald", cov_type = cov_type) |> print()
}

cov_types <- c("robust", "naive", "bias-corrected", "df-adjusted")
for(cov_type in cov_types) {
  anova(object2, test = "score", cov_type = cov_type) |> print()
  anova(object0, object2, test = "score", cov_type = cov_type) |> print()
  anova(object2, object0, object1, test = "score", cov_type = cov_type) |> print()
}

pmethods <- c("rao-scott", "satterthwaite")
for(pmethod in pmethods) {
  for(cov_type in cov_types) {
    anova(object2, test = "working-wald", cov_type = cov_type, pmethod = pmethod) |> print()
    anova(object0, object2, test = "working-wald", cov_type = cov_type, pmethod = pmethod) |> print()
    anova(object2, object0, object1, test = "working-wald", cov_type = cov_type, pmethod = pmethod) |> print()
  }
}


for(pmethod in pmethods) {
  for(cov_type in cov_types) {
    anova(object2, test = "working-score", cov_type = cov_type, pmethod = pmethod) |> print()
    anova(object0, object2, test = "working-score", cov_type = cov_type, pmethod = pmethod) |> print()
    anova(object2, object0, object1, test = "working-score", cov_type = cov_type, pmethod = pmethod) |> print()
  }
}


for(pmethod in pmethods) {
  for(cov_type in cov_types) {
    anova(object2, test = "working-lrt", cov_type = cov_type, pmethod = pmethod) |> print()
    anova(object0, object2, test = "working-lrt", cov_type = cov_type, pmethod = pmethod) |> print()
    anova(object2, object0, object1, test = "working-lrt", cov_type = cov_type, pmethod = pmethod) |> print()
  }
}




add1(object0, scope =  ~ period * treatment, test = "wald", cov_type = "robust")
add1(object0, scope =  ~ period * treatment, test = "working-wald")
add1(object0, scope =  ~ period * treatment, test = "score")
add1(object0, scope =  ~ period * treatment, test = "working-score")
add1(object0, scope =  ~ period * treatment, test = "working-lrt")


drop1(object0, test = "wald")
drop1(object0, test = "working-wald")
drop1(object0, test = "score")
drop1(object0, test = "working-score")
drop1(object0, test = "working-lrt")


drop1(object1, test = "wald")
drop1(object1, test = "working-wald")
drop1(object1, test = "score")
drop1(object1, test = "working-score")
drop1(object1, test = "working-lrt")


drop1(object2, test = "wald")
drop1(object2, test = "working-wald")
drop1(object2, test = "score")
drop1(object2, test = "working-score")
drop1(object2, test = "working-lrt")

drop1(object2, scope =  ~ 1, test = "wald")
drop1(object2, scope =  ~ 1, test = "working-wald")
drop1(object2, scope =  ~ 1, test = "score")
drop1(object2, scope =  ~ 1, test = "working-score")
drop1(object2, scope =  ~ 1, test = "working-lrt")


