fit_geewa_pois_indep <- geewa(
  formula = seizures ~ treatment + lnbaseline + lnage,
  data = test_data$epilepsy,
  id = id,
  family = poisson(link = "log"),
  corstr = "independence",
  method = "gee"
)

fit_geewa_pois_exch <- geewa(
  formula = seizures ~ treatment + lnbaseline + lnage,
  data = test_data$epilepsy,
  id = id,
  family = poisson(link = "log"),
  corstr = "exchangeable",
  method = "gee"
)

fit_geewa_bin_exch <- geewa_binary(
  formula = ecg ~ period + treatment,
  data = test_data$cerebrovascular,
  id = id,
  link = "logit",
  orstr = "exchangeable",
  method = "gee"
)
