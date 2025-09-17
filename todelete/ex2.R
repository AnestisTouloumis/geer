devtools::load_all()
N <- seq(10, 10000)
results1 <- rep(NA, length(N))
results2 <- rep(NA, length(N))
results3 <- rep(NA, length(N))

for(i in 1:length(N)) {
  x <- rep(c(0, 1), N[i])
  y <- rep(c(0, 1), N[i])
  id <- rep(1:N[i], each = 2)
  fit1 <-
    geewa(formula = y ~ x,
          id = id,
          corstr = "exchangeable",
          family = binomial(link = "logit"),
          method = "pgee-jeffreys",
          use_p = FALSE,
          phi_fixed = TRUE,
          control = geer_control(step_maxit = 1, tolerance = 1e-3))
  results1[i] <- fit1$iter
  results2[i] <- fit1$alpha
  results3[i] <- fit1$fitted[1]
  }
