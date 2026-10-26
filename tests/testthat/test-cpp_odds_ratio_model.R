list2env(get_cached_test_fixtures("cpp_odds_ratio_model"), environment())

test_that("Plackett covariance matches the stable closed form for moderate odds ratios", {
  for (mu in list(c(0.3, 0.6), c(0.05, 0.9), c(0.5, 0.5), c(0.2, 0.2))) {
    for (psi in c(0.02, 0.2, 0.7, 1.5, 4, 50, 1e3)) {
      V <- or_V(mu, psi)
      expect_rel(V[1, 2], ref_plackett(mu[1], mu[2], psi) - prod(mu),
        tolerance = 1e-8,
        label = sprintf("mu=(%g,%g) psi=%g", mu[1], mu[2], psi)
      )
      expect_equal(diag(V), mu * (1 - mu), tolerance = 1e-12)
    }
  }
})

test_that("weights enter the Plackett covariance as V / sqrt(w_i w_j)", {
  mu <- c(0.3, 0.6)
  psi <- 3
  V1 <- or_V(mu, psi)
  Vw <- or_V(mu, psi, w = c(2, 0.5))
  expect_equal(Vw, V1 / sqrt(outer(c(2, 0.5), c(2, 0.5))), tolerance = 1e-12)
})

test_that("Plackett covariance satisfies the Frechet bounds and approaches them", {
  for (mu in list(c(0.3, 0.6), c(0.7, 0.8), c(0.05, 0.9))) {
    b <- ref_frechet_cov(mu[1], mu[2])
    for (psi in c(1e-8, 1e-4, 0.5, 2, 1e4, 1e8)) {
      v12 <- or_V(mu, psi)[1, 2]
      expect_gte(v12, b[1] - 1e-12)
      expect_lte(v12, b[2] + 1e-12)
    }
    expect_equal(or_V(mu, 1e8)[1, 2], b[2], tolerance = 1e-6)
    expect_equal(or_V(mu, 1e-8)[1, 2], b[1], tolerance = 1e-6)
  }
})

test_that("Plackett covariance is symmetric in (p, q), monotone in psi and flips under complementing", {
  for (mu in list(c(0.3, 0.6), c(0.1, 0.45))) {
    for (psi in c(0.3, 2.5, 20)) {
      v <- or_V(mu, psi)[1, 2]
      expect_equal(or_V(rev(mu), psi)[1, 2], v, tolerance = 1e-12)
      expect_equal(or_V(c(mu[1], 1 - mu[2]), 1 / psi)[1, 2], -v, tolerance = 1e-9)
    }
  }
  grid <- sapply(c(0.1, 0.5, 1.2, 3, 10), function(psi) or_V(c(0.3, 0.6), psi)[1, 2])
  expect_true(all(diff(grid) > 0))
})

test_that("the implied joint table reproduces the working odds ratio", {
  mu <- c(0.35, 0.55)
  for (psi in c(0.1, 0.6, 1.7, 9)) {
    p11 <- or_V(mu, psi)[1, 2] + prod(mu)
    p10 <- mu[1] - p11
    p01 <- mu[2] - p11
    p00 <- 1 - mu[1] - mu[2] + p11
    expect_equal(p11 * p00 / (p10 * p01), psi, tolerance = 1e-9)
  }
})

test_that("Plackett covariance is relatively accurate for odds ratios close to one", {
  for (mu in list(c(0.3, 0.6), c(1e-3, 0.5))) {
    for (d in c(1e-4, 1e-5, 1e-6, 1e-7)) {
      psi <- 1 + d
      expect_rel(or_V(mu, psi)[1, 2], ref_plackett(mu[1], mu[2], psi) - prod(mu),
        tolerance = 1e-8,
        label = sprintf("mu=(%g,%g) psi-1=%g", mu[1], mu[2], d)
      )
    }
  }
})

test_that("fits are accurate and reach tight tolerance when the working odds ratios are near one", {
  d <- ref_sim("or", times8, link = "logit", alpha = rep(1 + 1e-7, 6), seed = 21)
  for (m in c("gee", "brgee-robust", "brgee-naive", "pgee-jeffreys")) {
    fit <- ref_fit(d, m, geer_control(tolerance = 1e-10, maxiter = 200))
    expect_true(fit$converged, label = m)
  }
})

test_that("get_marginalized_odds_ratios_or equals explicit weighted 2x2 tables with the Haldane constant", {
  withr::local_seed(2)
  times <- rep(list(1:4, c(1, 2, 4), 2:4, c(1, 3, 4)), 15)
  d <- ref_sim("or", times,
    link = "logit", alpha = rep(2, 6), weights = c(1, 2, 0.5, 1.5),
    seed = 61
  )
  ord <- order(d$id, d$tt)
  for (add in c(0.5, 0.1)) {
    for (str in c("unstructured", "exchangeable")) {
      got <- as.numeric(get_marginalized_odds_ratios_or(
        d$y[ord], d$id[ord], d$tt[ord], d$w[ord], add, str
      ))
      expected <- ref_or_counts(
        d$y[ord], d$id[ord], d$tt[ord], d$w[ord], add, 4,
        str == "exchangeable"
      )
      expect_equal(got, expected, tolerance = 1e-12, label = paste(add, str))
    }
  }
})

test_that("pairs of visits that are never observed together do not dilute the exchangeable odds ratio", {
  withr::local_seed(1)
  N <- 200
  u <- rbinom(N, 1, 0.5)
  pr <- ifelse(u == 1, 0.8, 0.2)
  y <- c(as.vector(t(cbind(rbinom(N, 1, pr), rbinom(N, 1, pr)))), rbinom(N, 1, 0.5))
  id <- c(rep(1:N, each = 2), N + 1:N)
  tt <- c(rep(1:2, N), rep(3L, N))
  o <- order(id, tt)
  un <- as.numeric(get_marginalized_odds_ratios_or(
    y[o], id[o], tt[o], rep(1, length(y)), 0.5, "unstructured"
  ))
  ex <- as.numeric(get_marginalized_odds_ratios_or(
    y[o], id[o], tt[o], rep(1, length(y)), 0.5, "exchangeable"
  ))
  expect_equal(un[2:3], c(1, 1))
  expect_equal(ex, rep(un[1], 3), tolerance = 1e-12)
})

test_that("bias reduction handles strong negative association", {
  respiratory <- test_data$respiratory
  fit <- geewa_binary(
    formula = status ~ 1,
    data = respiratory[respiratory$visit <= 2, ],
    id = id,
    repeated = visit,
    orstr = "fixed",
    alpha_vector = 0.05,
    method = "brgee-robust"
  )
  expect_true(fit$converged)
})

test_that("rare-event covariances retain Plackett association", {
  for (probability in c(1e-12, 1e-9, 1 - 1e-9)) {
    expected <- ref_plackett(
      min(probability, 1 - probability),
      min(probability, 1 - probability), 1e12
    ) - min(probability, 1 - probability)^2
    covariance <- or_V(rep(probability, 2L), 1e12)[1L, 2L]
    expect_equal(covariance / expected, 1, tolerance = 1e-6)
    expect_equal(
      compute_plackett_covariance(probability, probability, 1e12),
      covariance, tolerance = 1e-12
    )
  }
})

test_that("extreme finite odds ratios have finite limiting covariances", {
  for (odds_ratio in c(1e-200, 1e200)) {
    bounds <- ref_frechet_cov(0.3, 0.6)
    expected <- bounds[if (odds_ratio < 1) 1L else 2L]
    expect_equal(or_V(c(0.3, 0.6), odds_ratio)[1L, 2L], expected,
                 tolerance = 1e-12)
    expect_equal(compute_plackett_covariance(0.3, 0.6, odds_ratio),
                 expected, tolerance = 1e-12)
    root <- sqrt(odds_ratio)
    expect_equal(compute_plackett_covariance(0.5, 0.5, odds_ratio),
                 (root - 1) / (4 * (root + 1)), tolerance = 1e-12)
  }
  score <- as.numeric(get_estimating_equations_or(
    c(0, 1), diag(2L), c(1, 1), c(1, 2), c(1, 1),
    "identity", c(0.3, 0.6), c(0.3, 0.6), 1e200
  ))
  expect_equal(score, c(-10 / 3, 10 / 3), tolerance = 1e-6)
})

test_that("rare-event covariance derivatives match an independent Jacobian", {
  probabilities <- c(1e-9, 2e-9)
  beta <- qlogis(probabilities)
  response <- c(0, 1)
  quantities <- get_influence_quantities_or(
    response, diag(2L), c(1, 1), c(1, 2), c(1, 1), "logit",
    probabilities, beta, 1e12
  )
  expected <- -ref_jac(function(value) {
    compute_reference_or_score(value, response, 1e12)
  }, beta)
  expect_equal(quantities$sensitivity, expected, tolerance = 1e-6)
})

test_that("Plackett covariance preserves degenerate probability boundaries", {
  for (probability in c(0, 1)) {
    expect_identical(
      compute_plackett_covariance(probability, 0.4, 5), 0
    )
    expect_identical(
      compute_plackett_covariance(0.4, probability, 5), 0
    )
  }
})

test_that("near-complementary margins retain a valid Plackett covariance", {
  row <- 1 - 1e-16
  column <- 1e-16
  covariance <- compute_plackett_covariance(row, column, 1e-200)
  reversed <- compute_plackett_covariance(column, row, 1e-200)
  expect_equal(covariance / -9.9999999999999986e-17, 1,
               tolerance = 1e-12)
  expect_equal(covariance / reversed, 1, tolerance = 1e-12)
  expect_gte(covariance + row * column, 0)
  reference <- test_data$plackett_tail
  actual <- mapply(compute_plackett_covariance, reference$row,
                   reference$column, reference$odds_ratio)
  expect_lt(max(abs(actual / reference$covariance - 1)), 1e-12)
})

test_that("negative-association derivatives match an independent Jacobian", {
  probabilities <- c(0.2, 0.8)
  beta <- qlogis(probabilities)
  response <- c(0, 1)
  quantities <- get_influence_quantities_or(
    response, diag(2L), c(1, 1), c(1, 2), c(1, 1), "logit",
    probabilities, beta, 1e-4
  )
  expected <- -ref_jac(function(value) {
    compute_reference_or_score(value, response, 1e-4)
  }, beta, h = 1e-6)
  expect_equal(quantities$sensitivity, expected, tolerance = 1e-6)
})
