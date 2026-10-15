pearson_ref <- function(d, b) {
  mu <- d$lk$linkinv(drop(d$X %*% b))
  (d$y - mu) * sqrt(d$w / d$family$variance(mu))
}

alpha_ref <- function(d, r, phi, p, struct, m = 1L) {
  cl <- ref_clusters(d)
  n_occ <- d$tmax
  pairs <- do.call(rbind, lapply(cl, function(idx) {
    n <- length(idx)
    if (n < 2) {
      return(NULL)
    }
    jk <- t(combn(n, 2))
    data.frame(
      s = d$tt[idx][jk[, 1]], t = d$tt[idx][jk[, 2]],
      prod = r[idx][jk[, 1]] * r[idx][jk[, 2]]
    )
  }))
  pairs$lag <- pairs$t - pairs$s
  den <- function(n) (n - p) * phi
  switch(struct,
    exchangeable = sum(pairs$prod) / den(nrow(pairs)),
    ar1 = sum(pairs$prod[pairs$lag == 1]) / den(sum(pairs$lag == 1)),
    toeplitz = sapply(
      seq_len(n_occ - 1),
      function(l) sum(pairs$prod[pairs$lag == l]) / den(sum(pairs$lag == l))
    ),
    "m-dependent" = sapply(
      seq_len(m),
      function(l) sum(pairs$prod[pairs$lag == l]) / den(sum(pairs$lag == l))
    ),
    unstructured = unlist(lapply(1:(n_occ - 1), function(s) {
      sapply((s + 1):n_occ, function(t) {
        sel <- pairs$s == s & pairs$t == t
        sum(pairs$prod[sel]) / den(sum(sel))
      })
    }))
  )
}

times8 <- list(1:3, c(1, 2, 4), 2:4, 1:4, c(1, 4), 3:4, 1:4, 1:3)

test_that("get_pearson_residuals equals (y - mu) sqrt(w / V(mu)) for every variance function", {
  withr::local_seed(1)
  y <- runif(8, 0.1, 0.9)
  mu <- runif(8, 0.1, 0.9)
  w <- runif(8, 0.5, 3)
  for (fam in list(binomial(), poisson(), gaussian(), Gamma(), inverse.gaussian())) {
    expect_equal(drop(get_pearson_residuals(fam$family, y, mu, w)),
      (y - mu) * sqrt(w / fam$variance(mu)),
      tolerance = 1e-14, label = fam$family
    )
  }
})

test_that("estimated phi and alpha reproduce the moment formulas at the fitted mean", {
  ctl <- geer_control(tolerance = 1e-12, maxiter = 200)
  for (cf in list(
    list(poisson(), "log", 1L), list(binomial(), "logit", 1L),
    list(Gamma(), "log", 1L)
  )) {
    d <- ref_sim("cc", rep(times8, 10), cf[[1]], cf[[2]], "exchangeable", 0.3, 1, seed = 51)
    for (struct in c("exchangeable", "ar1", "toeplitz", "m-dependent", "unstructured")) {
      Mv <- if (struct == "m-dependent") 2L else 1L
      fit <- geewa(y ~ x1 + x2,
        data = d$df, id = id, repeated = visit, family = d$family,
        corstr = struct, Mv = Mv, control = ctl
      )
      b <- unname(coef(fit))
      r <- pearson_ref(d, b)
      p <- length(b)
      phi <- sum(r^2) / (length(r) - p)
      expect_equal(fit$phi, phi, tolerance = 1e-8, label = paste(cf[[1]]$family, struct, "phi"))
      expect_equal(as.vector(fit$alpha), alpha_ref(d, r, phi, p, struct, Mv),
        tolerance = 1e-8,
        label = paste(cf[[1]]$family, struct, "alpha")
      )
    }
  }
})

test_that("get_correlation_matrix builds the documented patterns", {
  gc <- get_correlation_matrix
  expect_equal(gc("exchangeable", 0.3, 4L), ref_corr("exchangeable", 0.3, 4), tolerance = 1e-15)
  expect_equal(gc("ar1", 0.5, 4L), ref_corr("ar1", 0.5, 4), tolerance = 1e-15)
  expect_equal(gc("independence", numeric(0), 3L), diag(3), tolerance = 0)
  a <- c(0.1, 0.2, 0.3, 0.15, 0.25, 0.05)
  expect_equal(gc("unstructured", a, 4L), ref_corr("unstructured", a, 4), tolerance = 1e-15)
  expect_equal(gc("fixed", a, 4L)[1, 4], 0.3)
  expect_equal(gc("fixed", a, 4L)[2, 4], 0.25)
  expect_equal(gc("toeplitz", c(0.5, 0.2, 0.1), 4L)[1, ], c(1, 0.5, 0.2, 0.1), tolerance = 1e-15)
  expect_equal(gc("m-dependent", c(0.5, 0.2), 4L)[1, ], c(1, 0.5, 0.2, 0), tolerance = 1e-15)
  expect_equal(gc("exchangeable", 0.3, 1L), matrix(1), tolerance = 0)
})

test_that("get_correlation_matrix rejects inadmissible or non-finite parameters", {
  gc <- get_correlation_matrix
  expect_error(gc("exchangeable", -0.6, 3L), "must lie in")
  expect_error(gc("exchangeable", 1, 3L), "must lie in")
  expect_error(gc("ar1", 1, 3L), "must lie in")
  expect_error(gc("ar1", -1.5, 3L), "must lie in")
  expect_error(gc("unstructured", c(0.9, -0.9, 0.9), 3L), "positive definite")
  expect_error(gc("exchangeable", NaN, 3L))
  expect_error(gc("ar1", NA_real_, 3L))
  expect_error(gc("exchangeable", numeric(0), 3L))
  expect_error(gc("ar1", numeric(0), 3L))
})

test_that("nuisance estimation errors are informative when association is not estimable", {
  df <- data.frame(id = 1:30, visit = 1L, x = seq(-1, 1, length.out = 30))
  df$y <- rpois(30, 3)
  expect_error(
    geewa(y ~ x, data = df, id = id, repeated = visit, family = poisson(), corstr = "exchangeable"),
    "cannot be estimated"
  )
  expect_error(
    geewa(y ~ x, data = df, id = id, repeated = visit, family = poisson(), corstr = "toeplitz"),
    "Toeplitz correlation cannot be estimated"
  )
})
