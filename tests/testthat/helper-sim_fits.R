influence_sim <- ref_sim(
  "cc",
  rep(list(1:4), 25),
  family = Gamma(),
  link = "log",
  weights = c(0.5, 2, 1, 1.5, 0.8)
)

fit_gamma_weighted <- geewa(
  formula = y ~ x1 + x2,
  data = influence_sim$df[-c(2, 11, 40), ],
  id = id,
  repeated = visit,
  family = Gamma(link = "log"),
  corstr = "unstructured",
  weights = w
)
