resolve_runs_order <- function(object, order_by, env = parent.frame()) {
  n_obs <- length(object$fitted.values)
  natural_order <- order(
    object$id,
    object$repeated,
    seq_len(n_obs)
  )

  if (is.character(order_by)) {
    if (length(order_by) != 1L || is.na(order_by) || !nzchar(order_by)) {
      stop("'order_by' must be a single non-empty character value or a numeric vector", call. = FALSE)
    }

    if (identical(order_by, "natural")) {
      return(list(
        index = natural_order,
        label = "natural cluster/repeated order",
        natural = TRUE
      ))
    }

    x_names <- colnames(object$x)
    if (identical(order_by, "fitted")) {
      key <- object$fitted.values
      label <- "fitted values"
    } else if (!is.null(x_names) && order_by %in% x_names) {
      key <- object$x[, order_by]
      label <- sprintf("model-matrix column '%s'", order_by)
    } else {
      key <- extract_runs_order_variable(object, order_by, env = env)
      label <- sprintf("variable '%s'", order_by)
    }
  } else if (is.numeric(order_by)) {
    if (length(order_by) != n_obs) {
      stop("numeric 'order_by' must have one value per fitted observation", call. = FALSE)
    }
    key <- as.numeric(order_by)
    label <- "supplied ordering vector"
  } else {
    stop("'order_by' must be a character value or a numeric vector", call. = FALSE)
  }

  if (anyNA(key) || any(!is.finite(key))) {
    stop("ordering values must all be finite and non-missing", call. = FALSE)
  }

  ordered_natural <- natural_order
  tie_break <- seq_along(ordered_natural)
  index <- ordered_natural[order(key[ordered_natural], tie_break)]

  list(index = index, label = label, natural = FALSE)
}


compute_runs_statistics <- function(residual_values,
                                    alternative = "two.sided") {
  alternative <- match.arg(alternative, c("two.sided", "less", "greater"))
  if (!is.numeric(residual_values) || length(residual_values) < 1L) {
    stop("residuals must be a non-empty numeric vector", call. = FALSE)
  }
  if (anyNA(residual_values) || any(!is.finite(residual_values))) {
    stop("residuals must all be finite and non-missing", call. = FALSE)
  }

  signs <- sign(residual_values)
  zero_no <- sum(signs == 0)
  signs <- signs[signs != 0]

  positive_no <- sum(signs > 0)
  negative_no <- sum(signs < 0)
  if (positive_no == 0L || negative_no == 0L) {
    stop(
      "the runs test requires at least one positive and one negative residual",
      call. = FALSE
    )
  }

  runs_no <- 1L
  if (length(signs) > 1L) {
    runs_no <- runs_no + sum(signs[-1L] != signs[-length(signs)])
  }

  n_nonzero <- positive_no + negative_no
  expected_runs <- 2 * positive_no * negative_no / n_nonzero + 1
  variance_runs <-
    2 * positive_no * negative_no *
    (2 * positive_no * negative_no - positive_no - negative_no) /
    (n_nonzero^2 * (n_nonzero - 1))

  if (!is.finite(variance_runs) || variance_runs <= 0) {
    stop(
      "the runs-test variance is not positive for this residual sign sequence",
      call. = FALSE
    )
  }

  statistic <- (runs_no - expected_runs) / sqrt(variance_runs)
  p_value <- switch(
    alternative,
    two.sided = 2 * stats::pnorm(-abs(statistic)),
    less = stats::pnorm(statistic),
    greater = stats::pnorm(statistic, lower.tail = FALSE)
  )

  list(
    runs = as.integer(runs_no),
    expected_runs = expected_runs,
    variance_runs = variance_runs,
    statistic = statistic,
    p_value = p_value,
    positive = as.integer(positive_no),
    negative = as.integer(negative_no),
    zero = as.integer(zero_no),
    nonzero = as.integer(n_nonzero)
  )
}
normalize_runs_nperm <- function(nperm) {
  if (is.null(nperm)) {
    return(NULL)
  }
  if (length(nperm) != 1L || !is.numeric(nperm) || is.na(nperm) ||
      !is.finite(nperm) || nperm < 1 || nperm != round(nperm)) {
    stop("'nperm' must be NULL or a single positive whole number", call. = FALSE)
  }
  as.integer(nperm)
}


compute_runs_permutation <- function(residual_values,
                                     cluster,
                                     runs,
                                     alternative,
                                     nperm) {
  signs <- sign(residual_values)
  keep <- signs != 0
  signs <- signs[keep]
  cluster <- cluster[keep]
  n_nonzero <- length(signs)
  if (n_nonzero != length(cluster)) {
    stop(
      "the runs test requires one cluster identifier per residual",
      call. = FALSE
    )
  }

  ## Signs are permuted within each cluster and the sequence is then read off
  ## in the original ordering. Sorting by cluster code and a fresh uniform key
  ## produces such a permutation in a single pass: both index vectors visit the
  ## clusters in the same order and with the same block sizes, so the
  ## assignment below only rearranges signs inside a cluster.
  cluster_codes <- match(cluster, unique(cluster))
  base_index <- order(cluster_codes)
  permuted_runs <- integer(nperm)
  for (i in seq_len(nperm)) {
    permuted_signs <- signs
    permuted_signs[base_index] <-
      signs[order(cluster_codes, stats::runif(n_nonzero))]
    permuted_runs[i] <- 1L +
      sum(permuted_signs[-1L] != permuted_signs[-n_nonzero])
  }

  ## The observed sequence is included in both tail counts, giving p-values
  ## that are valid Monte Carlo tests rather than plug-in estimates.
  p_less <- (1 + sum(permuted_runs <= runs)) / (nperm + 1)
  p_greater <- (1 + sum(permuted_runs >= runs)) / (nperm + 1)
  p_value <- switch(
    alternative,
    two.sided = min(1, 2 * min(p_less, p_greater)),
    less = p_less,
    greater = p_greater
  )

  degenerate <- length(unique(permuted_runs)) == 1L
  if (degenerate) {
    warning(
      "every within-cluster permutation gives the same number of runs, so the permutation test has no power here",
      call. = FALSE
    )
  }

  ## The permutation mean and standard deviation replace the exchangeable
  ## moments E(T) and V(T) when the statistic is standardized, so that the
  ## reported statistic and p-value refer to the same null distribution.
  list(
    p_value = p_value,
    nperm = nperm,
    degenerate = degenerate,
    mean = mean(permuted_runs),
    sd = stats::sd(permuted_runs)
  )
}
extract_runs_order_variable <- function(object, name, env = parent.frame()) {
  ## A variable that is not a model-matrix column may still be available in the
  ## data used to fit the model. It cannot be taken from object$data directly:
  ## rows are dropped by na.action and by subset, and the retained rows are then
  ## reordered by order(id, repeated) before being stored. The model frame is
  ## therefore rebuilt from the stored call with the variable appended to the
  ## formula, and the reconstruction is verified against the stored cluster and
  ## time indices before it is used.
  unknown_message <- sprintf(
    "unknown 'order_by' value '%s': use 'natural', 'fitted', a model-matrix column name, a variable in the data used to fit the model, or a numeric vector",
    name
  )
  if (is.null(object$call) || is.null(object$formula)) {
    stop(unknown_message, call. = FALSE)
  }
  extended_formula <- tryCatch(
    stats::update.formula(
      object$formula,
      stats::as.formula(paste0(". ~ . + `", name, "`"))
    ),
    error = function(e) NULL
  )
  if (is.null(extended_formula)) {
    stop(unknown_message, call. = FALSE)
  }
  model_call <- object$call
  model_call$formula <- extended_formula
  model_frame <- tryCatch(
    build_geer_model_frame(model_call, env = env),
    error = function(e) NULL
  )
  if (is.null(model_frame) || !(name %in% names(model_frame))) {
    stop(unknown_message, call. = FALSE)
  }
  if (nrow(model_frame) != object$obs_no) {
    stop(
      sprintf(
        "'order_by' variable '%s' cannot be aligned: including it changes the number of complete cases; supply an aligned numeric vector instead",
        name
      ),
      call. = FALSE
    )
  }
  id_repeated <- extract_geer_id_repeated(model_frame, nrow(model_frame))
  ord <- order(id_repeated$id, id_repeated$repeated)
  aligned <- !is.null(object$id) && !is.null(object$repeated) &&
    isTRUE(all.equal(
      as.numeric(id_repeated$id[ord]),
      as.numeric(object$id),
      check.attributes = FALSE
    )) &&
    isTRUE(all.equal(
      as.numeric(id_repeated$repeated[ord]),
      as.numeric(object$repeated),
      check.attributes = FALSE
    ))
  if (!aligned) {
    stop(
      sprintf(
        "'order_by' variable '%s' cannot be aligned with the observations stored in the fitted model; supply an aligned numeric vector instead",
        name
      ),
      call. = FALSE
    )
  }
  key <- model_frame[[name]]
  if (is.factor(key) || is.character(key) || is.logical(key)) {
    key <- as.numeric(factor(key))
  }
  if (!is.numeric(key)) {
    stop(
      sprintf("'order_by' variable '%s' must be numeric or a factor", name),
      call. = FALSE
    )
  }
  as.numeric(key)[ord]
}


compute_runs_exact_p <- function(positive_no,
                                 negative_no,
                                 runs,
                                 alternative) {
  ## Exact null distribution of the number of runs T in a random arrangement of
  ## n_p positive and n_n negative signs:
  ##   P(T = 2k)     = 2 C(n_p-1, k-1) C(n_n-1, k-1) / C(N, n_p)
  ##   P(T = 2k + 1) = [C(n_p-1, k) C(n_n-1, k-1) +
  ##                    C(n_p-1, k-1) C(n_n-1, k)] / C(N, n_p)
  ## Binomial coefficients are evaluated on the log scale, which keeps the
  ## denominator representable for large N; out-of-range terms give -Inf and
  ## hence contribute zero.
  n_nonzero <- positive_no + negative_no
  runs_values <- seq.int(2L, n_nonzero)
  half <- runs_values %/% 2L
  log_denominator <- lchoose(n_nonzero, positive_no)
  even <- runs_values %% 2L == 0L
  pmf <- numeric(length(runs_values))
  pmf[even] <- 2 * exp(
    lchoose(positive_no - 1L, half[even] - 1L) +
      lchoose(negative_no - 1L, half[even] - 1L) -
      log_denominator
  )
  pmf[!even] <-
    exp(
      lchoose(positive_no - 1L, half[!even]) +
        lchoose(negative_no - 1L, half[!even] - 1L) -
        log_denominator
    ) +
    exp(
      lchoose(positive_no - 1L, half[!even] - 1L) +
        lchoose(negative_no - 1L, half[!even]) -
        log_denominator
    )
  total <- sum(pmf)
  if (!is.finite(total) || abs(total - 1) > 1e-6) {
    stop(
      "the exact runs distribution could not be evaluated for these sign counts",
      call. = FALSE
    )
  }
  p_less <- sum(pmf[runs_values <= runs])
  p_greater <- sum(pmf[runs_values >= runs])
  p_value <- switch(
    alternative,
    two.sided = 2 * min(p_less, p_greater),
    less = p_less,
    greater = p_greater
  )
  min(1, p_value)
}


normalize_runs_exact <- function(exact, nperm) {
  if (!is.null(exact)) {
    if (length(exact) != 1L || !is.logical(exact) || is.na(exact)) {
      stop("'exact' must be NULL or a single non-missing logical value", call. = FALSE)
    }
    exact <- as.logical(exact)
    if (exact && !is.null(nperm)) {
      stop(
        "'exact' and 'nperm' must not both be requested: they select different reference distributions",
        call. = FALSE
      )
    }
  }
  exact
}
