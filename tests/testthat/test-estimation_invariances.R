all_methods <- c(
  "gee", "brgee-robust", "brgee-naive", "brgee-empirical", "bcgee-robust",
  "bcgee-naive", "bcgee-empirical", "pgee-jeffreys"
)

times6 <- list(1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3, 1:3)

test_that("estimators are equivariant under invertible linear reparametrization of the design", {
  A <- cbind(c(1, 0, 0), c(0, 1, 2), c(0, -1, 1))
  for (kind in c("cc", "or")) {
    d <- if (kind == "cc") {
      ref_sim("cc", times6, poisson(), "log", "exchangeable", 0.3, 1.3, seed = 31)
    } else {
      ref_sim("or", times6, link = "logit", alpha = c(1.5, 2.2, 0.9), seed = 31)
    }
    d2 <- d
    Z <- d$X %*% A
    d2$df$x1 <- Z[, 2]
    d2$df$x2 <- Z[, 3]
    for (m in all_methods) {
      f1 <- ref_fit(d, m, ctl_tight)
      f2 <- ref_fit(d2, m, ctl_tight)
      expect_equal(unname(coef(f2)), drop(solve(A, coef(f1))),
        tolerance = 1e-7,
        label = paste(kind, m)
      )
      expect_equal(unname(vcov(f2, "robust")), solve(A) %*% vcov(f1, "robust") %*% t(solve(A)),
        tolerance = 1e-7, label = paste(kind, m, "vcov")
      )
    }
  }
})

test_that("column order of the design does not matter", {
  d <- ref_sim("cc", times6, binomial(), "probit", "ar1", 0.4, 1, seed = 32)
  d2 <- d
  d2$df$x1 <- d$df$x2
  d2$df$x2 <- d$df$x1
  for (m in all_methods) {
    f1 <- ref_fit(d, m, ctl_tight)
    f2 <- ref_fit(d2, m, ctl_tight)
    expect_equal(unname(coef(f2))[c(1, 3, 2)], unname(coef(f1)), tolerance = 1e-7, label = m)
  }
})

test_that("flipping the binary response maps beta to -beta for symmetric links", {
  for (kind in c("cc", "or")) {
    for (link in c("logit", "probit")) {
      d <- if (kind == "cc") {
        ref_sim("cc", times6, binomial(), link, "exchangeable", 0.3, 1, seed = 33)
      } else {
        ref_sim("or", times6, link = link, alpha = c(1.5, 2.2, 0.9), seed = 33)
      }
      d2 <- d
      d2$df$y <- 1 - d$df$y
      for (m in all_methods) {
        f1 <- ref_fit(d, m, ctl_tight)
        f2 <- ref_fit(d2, m, ctl_tight)
        expect_equal(unname(coef(f2)), -unname(coef(f1)),
          tolerance = 1e-7,
          label = paste(kind, link, m)
        )
      }
    }
  }
})

test_that("results do not depend on cluster labels, cluster order or (for exchangeable) visit labels", {
  withr::local_seed(41)
  for (kind in c("cc", "or")) {
    d <- if (kind == "cc") {
      ref_sim("cc", times6, poisson(), "log", "exchangeable", 0.3, 1.3, seed = 34)
    } else {
      ref_sim("or", times6, link = "logit", alpha = rep(1.8, 3), seed = 34)
    }
    perm_ids <- sample(1000, length(unique(d$df$id)))
    d2 <- d
    d2$df$id <- perm_ids[d$df$id]
    d2$df <- d2$df[order(-d2$df$id), ]
    d3 <- d
    d3$df$visit <- c(3L, 1L, 2L)[d$df$visit]
    for (m in c("gee", "brgee-robust", "pgee-jeffreys")) {
      f1 <- ref_fit(d, m, ctl_tight)
      expect_equal(coef(ref_fit(d2, m, ctl_tight)), coef(f1),
        tolerance = 1e-9,
        label = paste(kind, m, "ids")
      )
      expect_equal(coef(ref_fit(d3, m, ctl_tight)), coef(f1),
        tolerance = 1e-9,
        label = paste(kind, m, "visits")
      )
    }
  }
})

test_that("duplicating every cluster leaves gee unchanged and halves the covariances", {
  for (kind in c("cc", "or")) {
    d <- if (kind == "cc") {
      ref_sim("cc", times6, poisson(), "log", "ar1", 0.4, 1.3, seed = 35)
    } else {
      ref_sim("or", times6, link = "logit", alpha = c(1.5, 2.2, 0.9), seed = 35)
    }
    d2 <- d
    dup <- d$df
    dup$id <- dup$id + max(dup$id)
    d2$df <- rbind(d$df, dup)
    f1 <- ref_fit(d, "gee", ctl_tight)
    f2 <- ref_fit(d2, "gee", ctl_tight)
    expect_equal(coef(f2), coef(f1), tolerance = 1e-9)
    expect_equal(vcov(f2, "naive"), vcov(f1, "naive") / 2, tolerance = 1e-9)
    expect_equal(vcov(f2, "robust"), vcov(f1, "robust") / 2, tolerance = 1e-9)
  }
})

test_that("only phi / w enters the correlation-based working covariance", {
  d <- ref_sim("cc", times6, poisson(), "log", "exchangeable", 0.3, 1.3,
    weights = c(1, 2, 0.5),
    seed = 36
  )
  d2 <- d
  d2$df$w <- 3 * d$df$w
  d2$phi <- 3 * d$phi
  for (m in all_methods) {
    f1 <- ref_fit(d, m, ctl_tight)
    f2 <- ref_fit(d2, m, ctl_tight)
    expect_equal(coef(f2), coef(f1), tolerance = 1e-8, label = m)
    expect_equal(vcov(f2, "robust"), vcov(f1, "robust"), tolerance = 1e-8, label = m)
  }
})

test_that("a Gaussian identity fit is equivariant in the response (scale and location)", {
  d <- ref_sim("cc", times6, gaussian(), "identity", "exchangeable", 0.3, 1, seed = 37)
  d2 <- d
  d2$df$y <- 5 * d$df$y + 2
  f1 <- geewa(y ~ x1 + x2,
    data = d$df, id = id, repeated = visit, corstr = "exchangeable",
    control = ctl_tight
  )
  f2 <- geewa(y ~ x1 + x2,
    data = d2$df, id = id, repeated = visit, corstr = "exchangeable",
    control = ctl_tight
  )
  expect_equal(unname(coef(f2)), unname(c(5 * coef(f1)[[1]] + 2, 5 * coef(f1)[-1])),
    tolerance = 1e-9
  )
  expect_equal(f2$phi, 25 * f1$phi, tolerance = 1e-9)
  expect_equal(f2$alpha, f1$alpha, tolerance = 1e-9)
})
