mcar_normal_initial_parameters <- function(x) {
  n <- nrow(x)
  p <- ncol(x)
  mu <- colMeans(x, na.rm = TRUE)
  x_imputed <- x
  for (j in seq_len(p)) {
    x_imputed[is.na(x_imputed[, j]), j] <- mu[j]
  }
  centered <- sweep(x_imputed, 2L, mu, FUN = "-")
  sigma <- crossprod(centered) / n

  for (j in seq_len(p)) {
    values <- x[!is.na(x[, j]), j]
    sigma[j, j] <- sum((values - mu[j])^2) / length(values)
  }

  sigma <- (sigma + t(sigma)) / 2
  scale <- max(c(diag(sigma), 1), na.rm = TRUE)
  ridge <- .Machine$double.eps^0.5 * scale
  attempts <- 0L
  while (inherits(try(chol(sigma), silent = TRUE), "try-error")) {
    sigma <- sigma + diag(ridge, p)
    ridge <- ridge * 10
    attempts <- attempts + 1L
    if (attempts > 12L) {
      stop(
        "initial covariance matrix for Little's MCAR test is singular; check for collinear or nearly constant variables",
        call. = FALSE
      )
    }
  }

  list(mu = mu, sigma = sigma)
}


mcar_normal_em <- function(x, maxit, tol) {
  n <- nrow(x)
  p <- ncol(x)
  initial <- mcar_normal_initial_parameters(x)
  mu <- initial$mu
  sigma <- initial$sigma
  missing <- is.na(x)
  pattern_key <- apply(missing, 1L, paste0, collapse = "")
  pattern_rows <- split(seq_len(n), pattern_key)

  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(maxit)) {
    sum_x <- numeric(p)
    sum_xx <- matrix(0, nrow = p, ncol = p)

    for (rows in pattern_rows) {
      observed <- which(!missing[rows[1L], ])
      missing_cols <- which(missing[rows[1L], ])
      n_pattern <- length(rows)

      if (!length(missing_cols)) {
        block <- x[rows, , drop = FALSE]
        sum_x <- sum_x + colSums(block)
        sum_xx <- sum_xx + crossprod(block)
        next
      }

      if (!length(observed)) {
        sum_x <- sum_x + n_pattern * mu
        sum_xx <- sum_xx + n_pattern * (sigma + tcrossprod(mu))
        next
      }

      sigma_oo <- sigma[observed, observed, drop = FALSE]
      sigma_mo <- sigma[missing_cols, observed, drop = FALSE]
      solved_oo <- tryCatch(
        solve(sigma_oo),
        error = function(e) {
          stop(
            "the observed-data covariance matrix became singular while fitting Little's MCAR test",
            call. = FALSE
          )
        }
      )
      regression <- sigma_mo %*% solved_oo
      conditional_cov <- sigma[missing_cols, missing_cols, drop = FALSE] -
        regression %*% sigma[observed, missing_cols, drop = FALSE]
      conditional_cov <- (conditional_cov + t(conditional_cov)) / 2

      observed_block <- x[rows, observed, drop = FALSE]
      centered_observed <- sweep(observed_block, 2L, mu[observed], FUN = "-")
      conditional_means <- sweep(
        centered_observed %*% t(regression),
        2L,
        mu[missing_cols],
        FUN = "+"
      )

      completed <- matrix(NA_real_, nrow = n_pattern, ncol = p)
      completed[, observed] <- observed_block
      completed[, missing_cols] <- conditional_means

      sum_x <- sum_x + colSums(completed)
      sum_xx <- sum_xx + crossprod(completed)
      sum_xx[missing_cols, missing_cols] <-
        sum_xx[missing_cols, missing_cols, drop = FALSE] +
        n_pattern * conditional_cov
    }

    mu_new <- sum_x / n
    sigma_new <- sum_xx / n - tcrossprod(mu_new)
    sigma_new <- (sigma_new + t(sigma_new)) / 2

    difference <- max(
      abs(mu_new - mu),
      abs(sigma_new - sigma)
    )
    parameter_scale <- max(1, abs(mu), abs(sigma))

    mu <- mu_new
    sigma <- sigma_new

    if (!all(is.finite(mu)) || !all(is.finite(sigma))) {
      stop("maximum-likelihood estimation for Little's MCAR test produced non-finite values", call. = FALSE)
    }

    if (difference <= tol * parameter_scale) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      sprintf("Little's MCAR test EM algorithm did not converge within %d iterations", maxit),
      call. = FALSE
    )
  }

  eigenvalues <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(1, max(abs(eigenvalues)))
  if (any(eigenvalues <= sqrt(.Machine$double.eps) * scale)) {
    stop(
      "the maximum-likelihood covariance matrix for Little's MCAR test is singular or nearly singular",
      call. = FALSE
    )
  }

  list(mu = mu, sigma = sigma, iterations = iter, converged = converged)
}


extract_geer_mcar_matrix <- function(object, data = NULL) {
  if (is.null(data)) {
    data <- object$data
  }
  if (is.null(data)) {
    stop(
      "the original data are not available in the fitted object; supply them through 'data'",
      call. = FALSE
    )
  }

  mcall <- object$call
  keep <- c("formula", "data", "subset", "id", "repeated", "na.action")
  mcall <- mcall[c(1L, match(keep, names(mcall), 0L))]
  mcall[[1L]] <- quote(stats::model.frame)
  mcall$formula <- stats::formula(object$terms)
  mcall$data <- quote(.geer_mcar_data)
  mcall$na.action <- quote(stats::na.pass)
  mcall$drop.unused.levels <- TRUE

  formula_env <- environment(stats::formula(object$terms))
  if (is.null(formula_env)) formula_env <- parent.frame()
  eval_env <- new.env(parent = formula_env)
  assign(".geer_mcar_data", data, envir = eval_env)
  model_frame <- tryCatch(
    eval(mcall, envir = eval_env, enclos = formula_env),
    error = function(e) {
      stop(
        sprintf("could not reconstruct the original data for Little's MCAR test: %s", conditionMessage(e)),
        call. = FALSE
      )
    }
  )

  response <- stats::model.response(model_frame)
  if (is.matrix(response) || is.data.frame(response)) {
    stop(
      "Little's MCAR test for a 'geer' object currently requires a univariate response at each repeated measurement",
      call. = FALSE
    )
  }
  if (!is.numeric(response)) {
    stop("the response used for Little's MCAR test must be numeric", call. = FALSE)
  }

  id_raw <- stats::model.extract(model_frame, "id")
  if (is.null(id_raw)) stop("'id' could not be recovered from the fitted model", call. = FALSE)
  if (anyNA(id_raw)) stop("'id' cannot contain missing values for Little's MCAR test", call. = FALSE)
  id <- as.numeric(factor(id_raw))

  repeated_raw <- stats::model.extract(model_frame, "repeated")
  if (is.null(repeated_raw)) {
    repeated <- stats::ave(id, id, FUN = seq_along)
    repeated_labels <- as.character(seq_len(max(repeated)))
  } else {
    if (anyNA(repeated_raw)) {
      stop("'repeated' cannot contain missing values for Little's MCAR test", call. = FALSE)
    }
    repeated_factor <- factor(repeated_raw)
    repeated <- as.numeric(repeated_factor)
    repeated_labels <- levels(repeated_factor)
  }

  if (length(response) != length(id) || length(response) != length(repeated)) {
    stop("response, 'id', and 'repeated' must have the same length", call. = FALSE)
  }
  if (any(duplicated(data.frame(id = id, repeated = repeated)))) {
    stop("'repeated' must identify unique measurements within each 'id'", call. = FALSE)
  }

  response_matrix <- matrix(
    NA_real_,
    nrow = max(id),
    ncol = max(repeated),
    dimnames = list(levels(factor(id_raw)), repeated_labels)
  )
  response_matrix[cbind(id, repeated)] <- as.numeric(response)
  response_matrix
}

