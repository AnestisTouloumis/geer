validate_mcar_homoscedasticity_matrix <- function(x) {
  if (is.data.frame(x)) {
    numeric_cols <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop(
        "all variables supplied to the Jamshidian-Jalal diagnostic must be numeric",
        call. = FALSE
      )
    }
    x <- as.matrix(x)
  } else if (is.matrix(x)) {
    if (!is.numeric(x)) {
      stop(
        "the matrix supplied to the Jamshidian-Jalal diagnostic must be numeric",
        call. = FALSE
      )
    }
  } else {
    stop(
      "'object' must be a fitted 'geer' object, numeric matrix, or numeric data frame",
      call. = FALSE
    )
  }

  storage.mode(x) <- "double"
  if (nrow(x) < 2L) {
    stop(
      "the Jamshidian-Jalal diagnostic requires at least two rows",
      call. = FALSE
    )
  }
  if (ncol(x) < 2L) {
    stop(
      "the Jamshidian-Jalal diagnostic requires at least two variables or repeated measurements",
      call. = FALSE
    )
  }
  if (any(is.infinite(x))) {
    stop(
      "data supplied to the Jamshidian-Jalal diagnostic cannot contain infinite values",
      call. = FALSE
    )
  }

  all_missing <- rowSums(!is.na(x)) == 0L
  if (any(all_missing)) {
    warning(
      sprintf(
        "%d row(s) with no observed values were omitted before the Jamshidian-Jalal diagnostic",
        sum(all_missing)
      ),
      call. = FALSE
    )
    x <- x[!all_missing, , drop = FALSE]
  }

  if (nrow(x) < 2L) {
    stop(
      "the Jamshidian-Jalal diagnostic requires at least two rows with observed data",
      call. = FALSE
    )
  }

  observed_no <- colSums(!is.na(x))
  labels <- colnames(x)
  if (is.null(labels)) labels <- as.character(seq_len(ncol(x)))
  if (any(observed_no < 2L)) {
    bad <- which(observed_no < 2L)
    stop(
      sprintf(
        paste0(
          "the Jamshidian-Jalal diagnostic requires at least two observed ",
          "values for variable(s) %s"
        ),
        paste(labels[bad], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  x
}

jj_missingness_key <- function(missing) {
  apply(
    missing,
    1L,
    function(z) paste0(as.integer(z), collapse = "")
  )
}

jj_pattern_information <- function(x, min_pattern_size) {
  missing <- is.na(x)
  key <- jj_missingness_key(missing)
  counts <- table(key)

  keep_keys <- names(counts)[counts >= min_pattern_size]
  omitted_keys <- setdiff(names(counts), keep_keys)
  keep <- key %in% keep_keys

  omitted <- NULL
  if (length(omitted_keys)) {
    omitted <- data.frame(
      pattern = omitted_keys,
      n = as.integer(counts[omitted_keys]),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    warning(
      sprintf(
        paste0(
          "%d missingness pattern(s) containing %d row(s) were omitted ",
          "because they had fewer than %d cases"
        ),
        length(omitted_keys),
        sum(omitted$n),
        min_pattern_size
      ),
      call. = FALSE
    )
  }

  x_used <- x[keep, , drop = FALSE]
  key_used <- key[keep]
  if (!nrow(x_used)) {
    stop(
      "no missingness patterns remain after applying 'min_pattern_size'",
      call. = FALSE
    )
  }

  pattern_levels <- unique(key_used)
  group <- match(key_used, pattern_levels)
  group_counts <- tabulate(group, nbins = length(pattern_levels))

  if (length(pattern_levels) < 2L) {
    stop(
      "the Jamshidian-Jalal diagnostic requires at least two retained missingness patterns",
      call. = FALSE
    )
  }

  pattern_matrix <- do.call(
    rbind,
    lapply(
      pattern_levels,
      function(z) as.integer(strsplit(z, "", fixed = TRUE)[[1L]])
    )
  )
  colnames(pattern_matrix) <- colnames(x_used)
  rownames(pattern_matrix) <- paste0("pattern", seq_along(pattern_levels))

  list(
    x = x_used,
    group = group,
    group_counts = group_counts,
    pattern_keys = pattern_levels,
    pattern_matrix = pattern_matrix,
    omitted_patterns = omitted
  )
}

jj_solve <- function(a, context) {
  out <- tryCatch(
    solve(a),
    error = function(e) NULL
  )
  if (is.null(out) || any(!is.finite(out))) {
    stop(
      sprintf("a covariance matrix was singular while %s", context),
      call. = FALSE
    )
  }
  out
}

jj_covariance_sqrt <- function(sigma, context) {
  sigma <- (sigma + t(sigma)) / 2
  eig <- eigen(sigma, symmetric = TRUE)
  scale <- max(1, max(abs(eig$values)))
  tolerance <- sqrt(.Machine$double.eps) * scale

  if (any(eig$values < -tolerance)) {
    stop(
      sprintf("a conditional covariance matrix was not positive semidefinite while %s", context),
      call. = FALSE
    )
  }

  values <- pmax(eig$values, 0)
  diag(sqrt(values), nrow = length(values)) %*% t(eig$vectors)
}

jj_normal_impute <- function(x, maxit, tol) {
  em <- mcar_normal_em(x, maxit = maxit, tol = tol)
  completed <- x
  missing <- is.na(x)
  pattern_key <- jj_missingness_key(missing)
  pattern_rows <- split(seq_len(nrow(x)), pattern_key)

  for (rows in pattern_rows) {
    missing_cols <- which(missing[rows[1L], ])
    if (!length(missing_cols)) next

    observed <- which(!missing[rows[1L], ])
    if (!length(observed)) {
      stop(
        "rows with no observed values cannot be imputed by the Jamshidian-Jalal diagnostic",
        call. = FALSE
      )
    }

    sigma_oo <- em$sigma[observed, observed, drop = FALSE]
    sigma_mo <- em$sigma[missing_cols, observed, drop = FALSE]
    regression <- sigma_mo %*% jj_solve(
      sigma_oo,
      "performing normal-theory imputation"
    )
    conditional_cov <- em$sigma[missing_cols, missing_cols, drop = FALSE] -
      regression %*% em$sigma[observed, missing_cols, drop = FALSE]
    conditional_cov <- (conditional_cov + t(conditional_cov)) / 2

    observed_block <- x[rows, observed, drop = FALSE]
    centered <- sweep(observed_block, 2L, em$mu[observed], FUN = "-")
    conditional_mean <- sweep(
      centered %*% t(regression),
      2L,
      em$mu[missing_cols],
      FUN = "+"
    )

    root <- jj_covariance_sqrt(
      conditional_cov,
      "performing normal-theory imputation"
    )
    innovations <- matrix(
      stats::rnorm(length(rows) * length(missing_cols)),
      nrow = length(rows),
      ncol = length(missing_cols)
    ) %*% root
    completed[rows, missing_cols] <- conditional_mean + innovations
  }

  list(
    data = completed,
    mu = em$mu,
    sigma = em$sigma,
    iterations = em$iterations,
    converged = em$converged
  )
}

jj_distribution_free_impute <- function(x) {
  complete <- stats::complete.cases(x)
  complete_data <- x[complete, , drop = FALSE]
  n_complete <- nrow(complete_data)
  p <- ncol(x)

  if (n_complete < 2L) {
    stop(
      "distribution-free imputation requires at least two complete cases",
      call. = FALSE
    )
  }

  mu <- colMeans(complete_data)
  sigma <- stats::cov(complete_data)
  eigenvalues <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(1, max(abs(eigenvalues)))
  if (any(eigenvalues <= sqrt(.Machine$double.eps) * scale)) {
    stop(
      paste0(
        "the complete-case covariance matrix is singular or nearly singular ",
        "for distribution-free imputation"
      ),
      call. = FALSE
    )
  }

  residuals <- sweep(complete_data, 2L, mu, FUN = "-") *
    sqrt(n_complete / (n_complete - 1))
  incomplete_rows <- which(!complete)
  sampled <- sample.int(
    n_complete,
    size = length(incomplete_rows),
    replace = TRUE
  )
  sampled_residuals <- residuals[sampled, , drop = FALSE]
  residual_position <- stats::setNames(
    seq_along(incomplete_rows),
    as.character(incomplete_rows)
  )

  completed <- x
  missing <- is.na(x)
  pattern_key <- jj_missingness_key(missing)
  pattern_rows <- split(incomplete_rows, pattern_key[incomplete_rows])

  for (rows in pattern_rows) {
    missing_cols <- which(missing[rows[1L], ])
    observed <- which(!missing[rows[1L], ])
    if (!length(missing_cols)) next
    if (!length(observed)) {
      stop(
        "rows with no observed values cannot be imputed by the Jamshidian-Jalal diagnostic",
        call. = FALSE
      )
    }

    s_oo <- sigma[observed, observed, drop = FALSE]
    a <- sigma[missing_cols, observed, drop = FALSE] %*%
      jj_solve(s_oo, "performing distribution-free imputation")

    observed_block <- x[rows, observed, drop = FALSE]
    predicted <- sweep(
      sweep(observed_block, 2L, mu[observed], FUN = "-") %*% t(a),
      2L,
      mu[missing_cols],
      FUN = "+"
    )

    residual_block <- sampled_residuals[
      residual_position[as.character(rows)],
      ,
      drop = FALSE
    ]
    innovation <- residual_block[, missing_cols, drop = FALSE] -
      residual_block[, observed, drop = FALSE] %*% t(a)

    completed[rows, missing_cols] <- predicted + innovation
  }

  list(
    data = completed,
    mu = mu,
    sigma = sigma,
    iterations = NA_integer_,
    converged = NA
  )
}

jj_impute <- function(x, imputation, maxit, tol) {
  p <- ncol(x)
  n_complete <- sum(stats::complete.cases(x))
  requested <- imputation

  if (identical(imputation, "distribution-free") &&
      (n_complete < 10L || n_complete < 2L * p)) {
    warning(
      sprintf(
        paste0(
          "distribution-free imputation requires at least 10 complete cases ",
          "and at least 2*p complete cases; only %d complete cases are available, ",
          "so normal-theory imputation is used"
        ),
        n_complete
      ),
      call. = FALSE
    )
    imputation <- "normal"
  }

  out <- if (identical(imputation, "distribution-free")) {
    jj_distribution_free_impute(x)
  } else {
    jj_normal_impute(x, maxit = maxit, tol = tol)
  }

  out$requested <- requested
  out$used <- imputation
  out$n_complete <- n_complete
  out
}

jj_legendre_values <- function(x) {
  z <- 2 * x - 1
  p0 <- rep(1, length(x))
  p1 <- z
  p2 <- (3 * z * p1 - p0) / 2
  p3 <- (5 * z * p2 - 2 * p1) / 3
  p4 <- (7 * z * p3 - 3 * p2) / 4

  cbind(
    sqrt(3) * p1,
    sqrt(5) * p2,
    sqrt(7) * p3,
    3 * p4
  )
}

jj_neyman_statistic <- function(x) {
  polynomials <- jj_legendre_values(x)
  sum(colSums(polynomials)^2) / length(x)
}

jj_neyman_p_value <- function(x, nrep, n_min) {
  statistic <- jj_neyman_statistic(x)
  n <- length(x)

  if (n >= n_min) {
    p_value <- stats::pchisq(statistic, df = 4, lower.tail = FALSE)
    return(list(statistic = statistic, p.value = p_value, simulated = FALSE))
  }

  simulated <- numeric(nrep)
  for (i in seq_len(nrep)) {
    simulated[i] <- jj_neyman_statistic(stats::runif(n))
  }
  p_value <- sum(simulated > statistic) / nrep
  if (p_value == 0) p_value <- 1 / nrep

  list(statistic = statistic, p.value = p_value, simulated = TRUE)
}

jj_hawkins_test <- function(
    completed,
    group,
    group_counts,
    nrep,
    n_min,
    test_uniformity = TRUE) {
  n <- nrow(completed)
  p <- ncol(completed)
  g <- length(group_counts)

  if (n - g - p <= 0L) {
    stop(
      "the Hawkins test requires n - number_of_patterns - number_of_variables > 0",
      call. = FALSE
    )
  }
  if (any(group_counts < 2L)) {
    stop("each retained missingness pattern must contain at least two cases", call. = FALSE)
  }

  pooled <- matrix(0, nrow = p, ncol = p)
  centered <- matrix(0, nrow = n, ncol = p)
  group_means <- matrix(NA_real_, nrow = g, ncol = p)

  for (i in seq_len(g)) {
    rows <- which(group == i)
    block <- completed[rows, , drop = FALSE]
    group_means[i, ] <- colMeans(block)
    centered[rows, ] <- sweep(block, 2L, group_means[i, ], FUN = "-")
    pooled <- pooled + (length(rows) - 1) * stats::cov(block)
  }
  pooled <- pooled / (n - g)
  pooled_inverse <- jj_solve(pooled, "calculating the Hawkins statistic")

  f_values <- numeric(n)
  uniform_values <- numeric(n)
  group_statistics <- numeric(g)
  group_p_values <- numeric(g)
  simulated <- logical(g)

  for (i in seq_len(g)) {
    rows <- which(group == i)
    ni <- group_counts[i]
    centered_block <- centered[rows, , drop = FALSE]
    v <- rowSums((centered_block %*% pooled_inverse) * centered_block)
    scaled_v <- ni * v
    denominator <- p * ((ni - 1) * (n - g) - scaled_v)

    if (any(denominator <= 0)) {
      tolerance <- sqrt(.Machine$double.eps) *
        max(1, max(abs(p * (ni - 1) * (n - g))))
      if (any(denominator < -tolerance)) {
        stop(
          paste0(
            "the Hawkins transformation produced a non-positive denominator; ",
            "the completed data may be too singular for this diagnostic"
          ),
          call. = FALSE
        )
      }
      denominator <- pmax(denominator, .Machine$double.eps)
    }

    f_i <- ((n - g - p) * scaled_v) / denominator
    a_i <- stats::pf(
      f_i,
      df1 = p,
      df2 = n - g - p,
      lower.tail = FALSE
    )
    f_values[rows] <- f_i
    uniform_values[rows] <- a_i

    if (test_uniformity) {
      neyman <- jj_neyman_p_value(a_i, nrep = nrep, n_min = n_min)
      group_statistics[i] <- neyman$statistic
      group_p_values[i] <- neyman$p.value
      simulated[i] <- neyman$simulated
    }
  }

  if (test_uniformity) {
    fisher <- -2 * sum(log(pmax(group_p_values, .Machine$double.xmin)))
    df <- 2L * g
    p_value <- stats::pchisq(fisher, df = df, lower.tail = FALSE)
  } else {
    fisher <- NA_real_
    df <- NA_integer_
    p_value <- NA_real_
    group_statistics[] <- NA_real_
    group_p_values[] <- NA_real_
    simulated[] <- FALSE
  }

  list(
    statistic = fisher,
    parameter = df,
    p.value = p_value,
    f.values = f_values,
    uniform.values = uniform_values,
    group.statistics = group_statistics,
    group.p.values = group_p_values,
    simulated = simulated,
    pooled.covariance = pooled,
    group.means = group_means
  )
}

jj_anderson_darling_test <- function(x, group, group_counts) {
  k <- length(group_counts)
  n <- length(x)
  if (k < 2L) {
    stop("the Anderson-Darling test requires at least two groups", call. = FALSE)
  }
  if (n <= 3L) {
    stop("the Anderson-Darling test requires more than three observations", call. = FALSE)
  }

  ordered_groups <- unlist(
    lapply(seq_len(k), function(i) which(group == i)),
    use.names = FALSE
  )
  x_grouped <- x[ordered_groups]
  cumulative_counts <- c(0L, cumsum(group_counts))

  x_sort <- sort(x_grouped)[seq_len(n - 1L)]
  distinct_index <- which(!duplicated(x_sort))
  counts <- c(distinct_index, length(x_sort) + 1L) - c(0L, distinct_index)
  h_j <- counts[seq.int(2L, length(distinct_index) + 1L)]
  h_n <- cumsum(h_j)
  z_j <- x_sort[distinct_index]

  ad_group <- numeric(k)
  for (i in seq_len(k)) {
    idx <- seq.int(cumulative_counts[i] + 1L, cumulative_counts[i + 1L])
    f_ij <- rowSums(outer(z_j, x_grouped[idx], FUN = "=="))
    m_ij <- cumsum(f_ij)
    numerator <- (n * m_ij - group_counts[i] * h_n)^2
    denominator <- h_n * (n - h_n)
    ad_group[i] <- sum(h_j * numerator / denominator) / group_counts[i]
  }

  statistic <- sum(ad_group) / n
  ad_group <- ad_group / n

  j_value <- sum(1 / group_counts)
  h_value <- sum(1 / seq_len(n - 1L))
  g_value <- 0
  if (n > 2L) {
    for (i in seq_len(n - 2L)) {
      g_value <- g_value +
        sum(1 / seq.int(i + 1L, n - 1L)) / (n - i)
    }
  }

  a <- (4 * g_value - 6) * (k - 1) + (10 - 6 * g_value) * j_value
  b <- (2 * g_value - 4) * k^2 + 8 * h_value * k +
    (2 * g_value - 14 * h_value - 4) * j_value -
    8 * h_value + 4 * g_value - 6
  c_value <- (6 * h_value + 2 * g_value - 2) * k^2 +
    (4 * h_value - 4 * g_value + 6) * k +
    (2 * h_value - 6) * j_value + 4 * h_value
  d <- (2 * h_value + 6) * k^2 - 4 * h_value * k

  variance <- (a * n^3 + b * n^2 + c_value * n + d) /
    ((n - 1) * (n - 2) * (n - 3))
  variance <- max(variance, 0)
  if (variance <= .Machine$double.eps) {
    stop(
      "the Anderson-Darling standardization variance is zero or numerically negligible",
      call. = FALSE
    )
  }

  standardized <- (statistic - (k - 1)) / sqrt(variance)

  alpha_grid <- c(0.25, 0.10, 0.05, 0.025, 0.01)
  b0 <- c(0.675, 1.281, 1.645, 1.960, 2.326)
  b1 <- c(-0.245, 0.250, 0.678, 1.149, 1.822)
  b2 <- c(-0.105, -0.305, -0.362, -0.391, -0.396)
  logit_alpha <- log((1 - alpha_grid) / alpha_grid)
  quantiles <- b0 + b1 / sqrt(k - 1) + b2 / (k - 1)
  interpolation_index <- if (standardized <= quantiles[3L]) 1:4 else 2:5
  spline_value <- stats::spline(
    quantiles[interpolation_index],
    logit_alpha[interpolation_index],
    xout = standardized
  )$y
  p_value <- 1 / (1 + exp(spline_value))

  list(
    statistic = statistic,
    standardized = standardized,
    variance = variance,
    p.value = p_value,
    group.statistics = ad_group
  )
}

jj_nonparametric_test <- function(hawkins, group, group_counts) {
  jj_anderson_darling_test(
    x = hawkins$f.values,
    group = group,
    group_counts = group_counts
  )
}

jj_with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))

  if (length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    stop("'seed' must be NULL or a single finite number", call. = FALSE)
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  force(code)
}

jj_interpret <- function(method, hawkins, nonparametric, alpha) {
  if (identical(method, "hawkins")) {
    if (hawkins$p.value <= alpha) {
      return(
        paste0(
          "The modified Hawkins test rejects the joint null of multivariate normality ",
          "and homogeneous covariance matrices. This alone does not distinguish ",
          "nonnormality from evidence against MCAR."
        )
      )
    }
    return(
      paste0(
        "The modified Hawkins test does not reject the joint null of multivariate ",
        "normality and homogeneous covariance matrices. There is no evidence against ",
        "MCAR from this screening diagnostic."
      )
    )
  }

  if (identical(method, "nonparametric")) {
    if (nonparametric$p.value <= alpha) {
      return(
        paste0(
          "The nonparametric Anderson-Darling test rejects homogeneity of the ",
          "pattern-specific covariance structure, providing evidence against MCAR ",
          "in the Jamshidian-Jalal screening framework."
        )
      )
    }
    return(
      paste0(
        "The nonparametric Anderson-Darling test does not reject homogeneity of the ",
        "pattern-specific covariance structure. There is no evidence against MCAR ",
        "from this screening diagnostic."
      )
    )
  }

  if (hawkins$p.value > alpha) {
    return(
      paste0(
        "The modified Hawkins test does not reject the joint null of multivariate ",
        "normality and homogeneous covariance matrices. There is no evidence against ",
        "MCAR from this screening diagnostic."
      )
    )
  }

  if (nonparametric$p.value > alpha) {
    return(
      paste0(
        "The modified Hawkins test rejects but the nonparametric Anderson-Darling ",
        "test does not. This is consistent with nonnormality rather than covariance ",
        "heterogeneity; there is no evidence against MCAR from the nonparametric screen."
      )
    )
  }

  paste0(
    "Both the modified Hawkins and nonparametric Anderson-Darling tests reject. ",
    "This provides evidence of covariance heterogeneity across missingness patterns ",
    "and therefore evidence against MCAR in the Jamshidian-Jalal screening framework."
  )
}
