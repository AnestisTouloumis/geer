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
  original_nrow <- nrow(x)
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

  row_index <- seq_len(original_nrow)
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
    row_index <- row_index[!all_missing]
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

  attr(x, "row.index") <- row_index
  attr(x, "original.nrow") <- original_nrow
  x
}

validate_mcar_homoscedasticity_imputed <- function(imputed_data, x) {
  if (is.data.frame(imputed_data)) {
    numeric_cols <- vapply(imputed_data, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop("all variables in 'imputed_data' must be numeric", call. = FALSE)
    }
    imputed_data <- as.matrix(imputed_data)
  } else if (!is.matrix(imputed_data) || !is.numeric(imputed_data)) {
    stop("'imputed_data' must be a numeric matrix or numeric data frame", call. = FALSE)
  }
  storage.mode(imputed_data) <- "double"

  original_nrow <- attr(x, "original.nrow")
  if (nrow(imputed_data) != original_nrow || ncol(imputed_data) != ncol(x)) {
    stop(
      sprintf(
        paste0(
          "'imputed_data' must have the same dimensions as the incomplete data: ",
          "%d rows and %d columns"
        ),
        original_nrow,
        ncol(x)
      ),
      call. = FALSE
    )
  }
  if (anyNA(imputed_data) || any(!is.finite(imputed_data))) {
    stop("'imputed_data' cannot contain missing or non-finite values", call. = FALSE)
  }

  imputed_data <- imputed_data[attr(x, "row.index"), , drop = FALSE]
  observed <- !is.na(x)
  discrepancy <- abs(imputed_data[observed] - x[observed])
  scale <- max(1, abs(x[observed]))
  if (any(discrepancy > sqrt(.Machine$double.eps) * scale)) {
    stop(
      "'imputed_data' must reproduce the observed values of the incomplete data",
      call. = FALSE
    )
  }
  dimnames(imputed_data) <- dimnames(x)
  imputed_data
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
    rows = which(keep),
    group = group,
    group_counts = group_counts,
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

jj_complete_case_moments <- function(x) {
  complete_data <- x[stats::complete.cases(x), , drop = FALSE]
  n_complete <- nrow(complete_data)

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

  list(mu = mu, sigma = sigma, residuals = residuals)
}

jj_imputation_setup <- function(x, imputation, maxit, tol) {
  p <- ncol(x)
  n_complete <- sum(stats::complete.cases(x))
  requested <- imputation

  # Distribution-free imputation requires max(10, 2 * p) complete cases, as in
  # the MissMech code and its warning message. Jamshidian, Jalal and Jansen
  # (2014, Sections 2.2 and 3.1) print min(10, 2 * p), which contradicts the
  # MissMech code and would admit as few as four complete cases when p = 2.
  if (identical(imputation, "distribution-free") &&
      (n_complete < 10L || n_complete < 2L * p)) {
    warning(
      sprintf(
        paste0(
          "distribution-free imputation requires at least 10 complete cases ",
          "and at least 2*p complete cases; only %d complete cases are available ",
          "among the retained rows, ",
          "so normal-theory imputation is used"
        ),
        n_complete
      ),
      call. = FALSE
    )
    imputation <- "normal"
  }

  if (identical(imputation, "distribution-free")) {
    moments <- jj_complete_case_moments(x)
    out <- list(
      mu = moments$mu,
      sigma = moments$sigma,
      residuals = moments$residuals,
      iterations = NA_integer_,
      converged = NA
    )
  } else {
    em <- mcar_normal_em(x, maxit = maxit, tol = tol)
    out <- list(
      mu = em$mu,
      sigma = em$sigma,
      residuals = NULL,
      iterations = em$iterations,
      converged = em$converged
    )
  }

  out$requested <- requested
  out$used <- imputation
  out$n_complete <- n_complete
  out
}

jj_normal_impute <- function(x, setup) {
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

    sigma_oo <- setup$sigma[observed, observed, drop = FALSE]
    sigma_mo <- setup$sigma[missing_cols, observed, drop = FALSE]
    regression <- sigma_mo %*% jj_solve(
      sigma_oo,
      "performing normal-theory imputation"
    )
    conditional_cov <- setup$sigma[missing_cols, missing_cols, drop = FALSE] -
      regression %*% setup$sigma[observed, missing_cols, drop = FALSE]
    conditional_cov <- (conditional_cov + t(conditional_cov)) / 2

    observed_block <- x[rows, observed, drop = FALSE]
    centered <- sweep(observed_block, 2L, setup$mu[observed], FUN = "-")
    conditional_mean <- sweep(
      centered %*% t(regression),
      2L,
      setup$mu[missing_cols],
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

  completed
}

jj_distribution_free_impute <- function(x, setup) {
  complete <- stats::complete.cases(x)
  n_complete <- nrow(setup$residuals)
  incomplete_rows <- which(!complete)
  sampled <- sample.int(
    n_complete,
    size = length(incomplete_rows),
    replace = TRUE
  )
  sampled_residuals <- setup$residuals[sampled, , drop = FALSE]
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

    s_oo <- setup$sigma[observed, observed, drop = FALSE]
    a <- setup$sigma[missing_cols, observed, drop = FALSE] %*%
      jj_solve(s_oo, "performing distribution-free imputation")

    observed_block <- x[rows, observed, drop = FALSE]
    predicted <- sweep(
      sweep(observed_block, 2L, setup$mu[observed], FUN = "-") %*% t(a),
      2L,
      setup$mu[missing_cols],
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

  completed
}

jj_draw_imputation <- function(x, setup) {
  completed <- if (identical(setup$used, "distribution-free")) {
    jj_distribution_free_impute(x, setup)
  } else {
    jj_normal_impute(x, setup)
  }
  if (anyNA(completed) || any(!is.finite(completed))) {
    stop("imputation produced missing or non-finite completed values", call. = FALSE)
  }
  completed
}

jj_neyman_statistics <- function(u) {
  z <- 2 * u - 1
  p1 <- z
  p2 <- (3 * z * p1 - 1) / 2
  p3 <- (5 * z * p2 - 2 * p1) / 3
  p4 <- (7 * z * p3 - 3 * p2) / 4

  (
    colSums(sqrt(3) * p1)^2 +
      colSums(sqrt(5) * p2)^2 +
      colSums(sqrt(7) * p3)^2 +
      colSums(3 * p4)^2
  ) / nrow(u)
}

jj_neyman_statistic <- function(x) {
  jj_neyman_statistics(matrix(x, ncol = 1L))
}

jj_neyman_null <- function(n, nrep, block_size = 10000L) {
  simulated <- numeric(nrep)
  done <- 0L
  while (done < nrep) {
    reps <- min(block_size, nrep - done)
    u <- matrix(stats::runif(n * reps), nrow = n, ncol = reps)
    simulated[done + seq_len(reps)] <- jj_neyman_statistics(u)
    done <- done + reps
  }
  simulated
}

jj_neyman_nulls <- function(group_counts, nrep, n_min) {
  lapply(
    group_counts,
    function(ni) if (ni < n_min) jj_neyman_null(ni, nrep) else NULL
  )
}

jj_neyman_p_value <- function(x, null) {
  statistic <- jj_neyman_statistic(x)

  if (is.null(null)) {
    p_value <- stats::pchisq(statistic, df = 4, lower.tail = FALSE)
    return(list(statistic = statistic, p.value = p_value, simulated = FALSE))
  }

  p_value <- (1 + sum(null >= statistic)) / (length(null) + 1)
  list(statistic = statistic, p.value = p_value, simulated = TRUE)
}

jj_hawkins_test <- function(
    completed,
    group,
    group_counts,
    neyman_nulls,
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
  if (test_uniformity && (!is.list(neyman_nulls) || length(neyman_nulls) != g)) {
    stop(
      "'neyman_nulls' must be a list with one element per retained missingness pattern",
      call. = FALSE
    )
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

    # (ni - 1)(n - g) - ni * v = (ni - 1)(n - g)(1 - u'W^{-1}u), where W is the
    # pooled within-pattern SSCP matrix and W - uu' is the SSCP matrix with the
    # case deleted, so the complement is nonnegative and vanishes only when
    # that deletion leaves a singular pooled covariance matrix.
    complement <- 1 - scaled_v / ((ni - 1) * (n - g))
    if (any(complement <= sqrt(.Machine$double.eps))) {
      stop(
        paste0(
          "the Hawkins transformation is undefined because deleting a case ",
          "leaves a singular pooled covariance matrix"
        ),
        call. = FALSE
      )
    }
    denominator <- p * ((ni - 1) * (n - g) - scaled_v)

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
      neyman <- jj_neyman_p_value(a_i, null = neyman_nulls[[i]])
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

# Upper-tail reference quantiles of the standardized k-sample Anderson-Darling
# statistic A2_kN (Scholz and Stephens, 1987, first version, not adjusted for
# ties), indexed by m = k - 1. Values are taken from ad.pval() in the kSamples
# package (Scholz and Zhu, version 1.2-12, GPL (>= 2)), where they were
# obtained by simulation with 2e6 replications and sample sizes of 500 per
# group. Rows correspond to jj_ad_probabilities and columns to jj_ad_m_grid.
jj_ad_m_grid <- c(1, 2, 3, 4, 6, 8, 10, Inf)

jj_ad_probabilities <- c(
  0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.025, 0.05,
  0.075, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.925, 0.95,
  0.975, 0.99, 0.9925, 0.995, 0.9975, 0.999, 0.99925, 0.9995, 0.99975,
  0.9999, 0.999925, 0.99995, 0.999975, 0.99999
)

jj_ad_quantiles <- matrix(
  c(
     -1.1954,  -1.5806,  -1.8172,  -2.0032,  -2.2526,  -2.4204,  -2.5283,  -4.2649,
     -1.1786,  -1.5394,  -1.7728,  -1.9426,  -2.1685,  -2.3288,  -2.4374,  -3.8906,
     -1.1660,  -1.5193,  -1.7462,  -1.9067,  -2.1260,  -2.2818,  -2.3926,  -3.7190,
     -1.1407,  -1.4659,  -1.6710,  -1.8105,  -2.0048,  -2.1356,  -2.2348,  -3.2905,
     -1.1253,  -1.4371,  -1.6314,  -1.7619,  -1.9396,  -2.0637,  -2.1521,  -3.0902,
     -1.0777,  -1.3503,  -1.5102,  -1.6177,  -1.7610,  -1.8537,  -1.9178,  -2.5758,
     -1.0489,  -1.2984,  -1.4415,  -1.5355,  -1.6625,  -1.7380,  -1.7936,  -2.3263,
     -0.9978,  -1.2098,  -1.3251,  -1.4007,  -1.4977,  -1.5555,  -1.5941,  -1.9600,
     -0.9417,  -1.1187,  -1.2090,  -1.2671,  -1.3382,  -1.3790,  -1.4050,  -1.6449,
     -0.8981,  -1.0491,  -1.1235,  -1.1692,  -1.2249,  -1.2552,  -1.2755,  -1.4395,
     -0.8598,  -0.9904,  -1.0513,  -1.0879,  -1.1317,  -1.1550,  -1.1694,  -1.2816,
     -0.7258,  -0.7938,  -0.8188,  -0.8312,  -0.8435,  -0.8471,  -0.8496,  -0.8416,
     -0.5966,  -0.6170,  -0.6177,  -0.6139,  -0.6073,  -0.5987,  -0.5941,  -0.5244,
     -0.4572,  -0.4383,  -0.4190,  -0.4033,  -0.3834,  -0.3676,  -0.3587,  -0.2533,
     -0.2966,  -0.2428,  -0.2078,  -0.1844,  -0.1548,  -0.1346,  -0.1224,   0.0000,
     -0.1009,  -0.0169,   0.0304,   0.0596,   0.0933,   0.1156,   0.1294,   0.2533,
      0.1571,   0.2635,   0.3169,   0.3480,   0.3823,   0.4038,   0.4166,   0.5244,
      0.5357,   0.6496,   0.6992,   0.7246,   0.7528,   0.7683,   0.7771,   0.8416,
      1.2255,   1.2989,   1.3202,   1.3254,   1.3305,   1.3286,   1.3257,   1.2816,
      1.5262,   1.5677,   1.5709,   1.5663,   1.5561,   1.5449,   1.5356,   1.4395,
      1.9633,   1.9430,   1.9190,   1.8975,   1.8641,   1.8389,   1.8212,   1.6449,
      2.7314,   2.5899,   2.5000,   2.4451,   2.3664,   2.3155,   2.2823,   1.9600,
      3.7825,   3.4425,   3.2582,   3.1423,   3.0036,   2.9101,   2.8579,   2.3263,
      4.1241,   3.7160,   3.4984,   3.3651,   3.2003,   3.0928,   3.0311,   2.4324,
      4.6044,   4.0847,   3.8348,   3.6714,   3.4721,   3.3453,   3.2777,   2.5758,
      5.4090,   4.7223,   4.4022,   4.1791,   3.9357,   3.7809,   3.6963,   2.8070,
      6.4954,   5.5823,   5.1456,   4.8657,   4.5506,   4.3275,   4.2228,   3.0902,
      6.8279,   5.8282,   5.3658,   5.0749,   4.7318,   4.4923,   4.3642,   3.1747,
      7.2755,   6.1970,   5.6715,   5.3642,   4.9991,   4.7135,   4.5945,   3.2905,
      8.1885,   6.8537,   6.2077,   5.8499,   5.4246,   5.1137,   4.9555,   3.4808,
      9.3061,   7.6592,   6.8500,   6.4806,   5.9919,   5.6122,   5.5136,   3.7190,
      9.6132,   7.9234,   7.1025,   6.6731,   6.1549,   5.8217,   5.7345,   3.7911,
     10.0989,   8.2395,   7.4326,   6.9567,   6.3908,   6.0110,   5.9566,   3.8906,
     10.8825,   8.8994,   7.8934,   7.4501,   6.9009,   6.4538,   6.2705,   4.0556,
     11.8537,   9.5482,   8.5568,   8.0283,   7.4418,   6.9524,   6.6195,   4.2649
  ),
  nrow = 35L,
  ncol = 8L,
  byrow = TRUE
)

# Follows kSamples::ad.pval(): the quantiles for each probability are
# interpolated to 1 / sqrt(m) by smoothing splines (spar = 0.4), and the upper
# log-odds are then fitted against the interpolated quantiles by a smoothing
# spline (spar = 0.25), with linear extrapolation beyond the tabulated range.
jj_ad_p_value <- function(standardized, m) {
  grid <- 1 / sqrt(jj_ad_m_grid)
  target <- 1 / sqrt(m)
  quantiles <- vapply(
    seq_len(nrow(jj_ad_quantiles)),
    function(i) {
      fit <- stats::smooth.spline(grid, jj_ad_quantiles[i, ], spar = 0.4)
      stats::predict(fit, target)$y
    },
    numeric(1)
  )
  upper <- 1 - jj_ad_probabilities
  fit <- stats::smooth.spline(quantiles, stats::qlogis(upper), spar = 0.25)
  log_odds <- stats::predict(fit, standardized)$y

  list(
    p.value = stats::plogis(log_odds),
    extrapolated = standardized < min(quantiles) || standardized > max(quantiles)
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
  reference <- jj_ad_p_value(standardized, m = k - 1)

  list(
    statistic = statistic,
    standardized = standardized,
    variance = variance,
    p.value = reference$p.value,
    extrapolated = reference$extrapolated,
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
