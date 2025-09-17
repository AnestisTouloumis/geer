devtools::load_all()
N <- 1000
x <- rep(c(0, 1), N)
y <- rep(c(0, 1), N)
id <- rep(1:N, each = 2)
or <- (1/2)/(N + 0.5)
p1 <- (1/2) / (N + 1)
p2 <- (N + 1/2)/(N + 1)

fit1 <-
  geewa(formula = y ~ x,
        id = id,
        corstr = "exchangeable",
        family = binomial(link = "logit"),
        method = "pgee-jeffreys",
        use_p = FALSE,
        phi_fixed = TRUE,
        control = geer_control(step_maxit = 1))

fitted(fit1)[1]

results1 <- matrix(NA, nrow = length(10:1000), 6)
results2 <- matrix(NA, nrow = length(10:1000), 6)
i = 1
for(N in 10:1000) {
  p1 <- (1/2) / (N + 1)
  p2 <- (N + 1/2)/(N + 1)
  x <- rep(c(0, 1), N)
  y <- rep(c(0, 1), N)
  id <- rep(1:N, each = 2)
  fit1 <-
    geewa(formula = y ~ x,
          id = id,
          corstr = "exchangeable",
          family = binomial(link = "logit"),
          method = "pgee-jeffreys",
          use_p = FALSE,
          phi_fixed = TRUE,
          beta_start = qlogis(p1/(1-p1)) * c(1, -2),
          control = geer_control(step_maxit = 1))
p0 <- (1/2) * (1/(N+1))
p1 <- (N + 1) / (2 * N^2 + 3 * N + 2)
p2 <- (2 * N^2 + 3 * N + 2)/(4 * (N^3 + 2 * N^2 + 2 * N + 1))
p3 <- (2 * N^3 + 4 * N^2 + 4 * N + 2)/(4 * N^4 + 10 * N^3 + 13 * N^2 + 10 * N + 4)
p4 <- (4 * N^4 + 10 * N^3 + 13 * N^2 + 10 * N + 4)/ (2 * (4 * N^5 + 12 * N^4 + 19 * N^3 + 19 * N^2 + 12 * N +4))
results1[i, ] <- c(p0, p1, p2, p3, p4, fit1$fitted[1])
results2[i, ] <- -c(p0/(1-p0), p1/(1 - p1), p2/(1 - p2), p3/(1 - p3), p4/(1 - p4), -fit1$alpha)
i <- i + 1
}

fit2 <-
  geewa(formula = y ~ x,
        id = id,
        corstr = "exchangeable",
        family = binomial(link = "logit"),
        method = "pgee-jeffreys",
        phi_fixed = FALSE,
        use_p = FALSE)

fit3 <-
  geewa_binary(formula = y ~ x,
               id = id,
               orstr = "exchangeable",
               link = "logit",
               method = "pgee-jeffreys")

fit4 <-
  geewa_binary(formula = y ~ x,
               id = id,
               orstr = "exchangeable",
               link = "logit",
               method = "pgee-jeffreys",
               control = geer_control(maxiter = 1, step_maxit = 1))

library("geessbin")

fit12 <- geessbin(formula = y ~ x,
                  id = id,
                  corstr = "exchangeable",
                  beta.method = "PGEE",
                  SE.method = "MB",
                  scale.fix = FALSE)
fit22 <- geessbin(formula = y ~ x,
                  id = id,
                  corstr = "exchangeable",
                  beta.method = "PGEE",
                  SE.method = "MB",
                  scale.fix = TRUE)

fit23 <- geessbin(formula = y ~ x,
                  id = id,
                  corstr = "independence",
                  beta.method = "PGEE",
                  SE.method = "MB",
                  scale.fix = FALSE)
