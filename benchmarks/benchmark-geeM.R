# =============================================================================
# Benchmark: geeM (geem) versus geer::geewa
# =============================================================================

library(geeM)
library(testthat)
bench_dir <- if (file.exists("benchmark-helpers.R")) "." else "benchmarks"
devtools::load_all(file.path(bench_dir, ".."))
source(file.path(bench_dir, "benchmark-helpers.R"))

geem_tol <- 1e-10

as_ref_geem <- function(ref) {
  s <- summary(ref)
  list(
    coefficients = coef(ref),
    naive = as.matrix(ref$naiv.var),
    robust = as.matrix(ref$var),
    fitted = fitted(ref),
    summary = list(cov_type = "robust", se = s$se.robust,
                   statistic = s$wald.test)
  )
}

fixed_from_geem <- function(args, ref) {
  geewa_fixed(args, as.matrix(ref$biggest.R.alpha), ref$phi)
}

check_geem <- function(label, ref, fit, offset = 0) {
  check_fit(label, as_ref_geem(ref), fit)
  check_fit(paste(label, "naive standard errors"),
            list(summary = list(cov_type = "naive",
                                se = summary(ref)$se.model)), fit)
  describe(paste(label, "relationships"), {
    it("summary() p-values = round(geer p-values, 8)", {
      p_geer <- summary(fit, cov_type = "robust")$coefficients[, "Pr(>|z|)"]
      expect_agree(summary(ref)$p, round(p_geer, 8))
    })
    it("predict() = geer linear predictor - offset", {
      expect_agree(in_geer_order(predict(ref), fit),
                   predict(fit, type = "link") - offset)
    })
  })
}


# =============================================================================
# Example 1: simulated negative binomial counts, poisson(log), ar1
# =============================================================================

generatedata <- function(beta, alpha, gamma, X, T, n) {
  mean.vec <- exp(crossprod(t(X), beta))
  y <- matrix(0, nrow = n, ncol = T)
  y[, 1] <- rnbinom(n, mu = mean.vec[1], size = mean.vec[1] / gamma)
  for (i in 1:n) {
    for (t in 2:T) {
      innovation.mean <- mean.vec[t] - alpha * (sqrt(mean.vec[t] *
                                                       mean.vec[t - 1]))
      I <- rnbinom(1, mu = innovation.mean, size = innovation.mean / gamma)
      first.shape <- alpha * sqrt(mean.vec[t] * mean.vec[t - 1]) / gamma
      second.shape <- mean.vec[t - 1] / gamma - first.shape
      u <- rbeta(1, shape1 = first.shape, shape2 = second.shape)
      a <- rbinom(1, size = y[i, t - 1], prob = u)
      y[i, t] <- a + I
    }
  }
  longform <- c(t(y))
  data.frame(count = longform, time = rep(X[, 2], times = n),
             subject = rep(c(1:n), each = T))
}
set.seed(1)
X <- cbind(rep(1, 5), c(-.5, -.25, 0, .25, .5))
testdat <- generatedata(beta = c(1, .5), alpha = .2, gamma = .5, X = X,
                        T = 5, n = 3000)

ref1 <- geem(count ~ time, id = subject, data = testdat, family = poisson,
             corstr = "ar1", tol = geem_tol)
fit1 <- fixed_from_geem(list(formula = count ~ time, id = quote(subject),
                             data = testdat, family = poisson(link = "log")),
                        ref1)
check_geem("Ex1 simulated counts poisson/log ar1", ref1, fit1)


# =============================================================================
# Example 2: ohio (geepack), binomial(logit), m-dependent with Mv = 1
# =============================================================================

data("ohio", package = "geepack")
args_ohio <- list(formula = resp ~ age + smoke + age:smoke, id = quote(id),
                  data = ohio)

ref2 <- geem(resp ~ age + smoke + age:smoke, id = id, data = ohio,
             family = binomial, corstr = "m-dep", Mv = 1, tol = geem_tol)
fit2 <- fixed_from_geem(c(args_ohio,
                          family = list(binomial(link = "logit"))), ref2)
check_geem("Ex2 ohio binomial/logit m-dependent(1)", ref2, fit2)


# =============================================================================
# Example 3: ohio (geepack), user-supplied cauchit link, m-dependent Mv = 1
# =============================================================================

LinkFun <- function(arg) qcauchy(arg)
InvLink <- function(arg) pcauchy(arg)
InvLinkDeriv <- function(arg) dcauchy(arg)
VarFun <- function(arg) arg * (1 - arg)
FunList <- list(LinkFun, VarFun, InvLink, InvLinkDeriv)

ref3 <- geem(resp ~ age + smoke + age:smoke, id = id, data = ohio,
             family = FunList, corstr = "m-dep", Mv = 1, tol = geem_tol)
fit3 <- fixed_from_geem(c(args_ohio,
                          family = list(binomial(link = "cauchit"))), ref3)
check_geem("Ex3 ohio binomial/cauchit m-dependent(1)", ref3, fit3)


# =============================================================================
# Example 4: seizure (geepack), poisson(log), exchangeable, with offset
# =============================================================================

data("seizure", package = "geepack")
seiz.l <- reshape(seizure,
                  varying = list(c("base", "y1", "y2", "y3", "y4")),
                  v.names = "y", times = 0:4, direction = "long")
seiz.l <- seiz.l[order(seiz.l$id, seiz.l$time), ]
seiz.l$t <- ifelse(seiz.l$time == 0, 8, 2)
seiz.l$x <- ifelse(seiz.l$time == 0, 0, 1)

ref4 <- geem(y ~ x + trt + x:trt + offset(log(t)), id = id, data = seiz.l,
             family = poisson, corstr = "exchangeable", tol = geem_tol)
fit4 <- fixed_from_geem(list(formula = y ~ x + trt + x:trt + offset(log(t)),
                             id = quote(id), data = seiz.l,
                             family = poisson(link = "log")), ref4)
check_geem("Ex4 seizure poisson/log exchangeable with offset", ref4, fit4,
           offset = in_geer_order(log(seiz.l$t), fit4))
