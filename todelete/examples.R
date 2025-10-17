devtools::load_all()


## Example 1: cerebrovascular
data("cerebrovascular")
fitted_model1 <-
  geewa_binary(formula = ecg ~ factor(treatment) * factor(period), id = id,
               repeated = period, link = "logit", data = cerebrovascular,
               orstr = "exchangeable", method = "gee")

## Example 2: cerebrovascular
data("cholecystectomy")
fitted_model2 <-
  geewa_binary(formula =
                 pain ~ I(treatment == "active") + age + I(gender=="female") + I(time>=5),
               id = id, repeated = time, link = "logit", data = cholecystectomy,
               orstr = "unstructured", method = "brgee-robust")

## Example 3: depression
data("depression")
fitted_model3 <-
  geewa(formula = score ~ treatment + visit + baseline, id = id, repeated = visit,
        family = gaussian(link = "identity"), data = depression, corstr = "exchangeable",
        method = "gee")

## Example 4: epilepsy
data("epilepsy")
epilepsy$visit4 <- with(epilepsy, I(visit == 4))
fitted_model4 <-
  geewa(formula = seizures ~ treatment*lnbaseline + lnage + visit4, id = id,
        repeated = visit, family = quasi(link = "log", variance = "mu^2"),
        data = epilepsy, corstr = "exchangeable", method = "brgee-robust")

## Example 5: leprosy
data("leprosy")
leprosy$periodn <- ifelse(leprosy$period == "pre", 0, 1)
leprosy$treatn <- ifelse(leprosy$treatment == "C", 0, 1)
fitted_model5 <-
  geewa(formula = bacilli ~ periodn + periodn:treatn, id = id, repeated = period,
        family = poisson(link = "log"), data = leprosy, corstr = "exchangeable",
        method = "gee")

## Example 5: respiratory
data("respiratory")
leprosy$periodn <- ifelse(leprosy$period == "pre", 0, 1)
leprosy$treatn <- ifelse(leprosy$treatment == "C", 0, 1)
fitted_model5 <-
  geewa(formula = bacilli ~ periodn + periodn:treatn, id = id, repeated = period,
        family = poisson(link = "log"), data = leprosy, corstr = "exchangeable",
        method = "gee")

## Example 6: respiratory
data("respiratory")
respiratory2 <-
  respiratory[respiratory$center == "C2", ]
fitted_model6 <-
  geewa_binary(formula =
                 status ~ I(treatment == "active")*gender + visit*age + baseline,
               id = id, repeated = visit, link = "probit", data = respiratory2,
               orstr = "unstructured", method = "pgee-jeffreys")

## Example 7: rinse
data("rinse")
fitted_model7 <-
    geewa(formula =
            score ~ (I(treatment == "A") + I(treatment == "B")) + baseline + factor(time),
      id = id, repeated = time, family = quasi(link = "log", variance = "mu^2"),
      data = rinse, corstr = "exchangeable", method = "gee")
