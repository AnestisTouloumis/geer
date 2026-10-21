compute_mcar_em_start <- function(x) {
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

  if (anyNA(mu) || !all(is.finite(sigma))) {
    stop(
      paste0(
        "initial covariance matrix for normal-theory EM estimation is ",
        "singular; check for collinear or nearly constant variables"
      ),
      call. = FALSE
    )
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
        paste0(
          "initial covariance matrix for normal-theory EM estimation is ",
          "singular; check for collinear or nearly constant variables"
        ),
        call. = FALSE
      )
    }
  }

  list(mu = mu, sigma = sigma)
}


fit_mcar_normal_em <- function(x, maxiter, tolerance) {
  n <- nrow(x)
  p <- ncol(x)
  initial <- compute_mcar_em_start(x)
  mu <- initial$mu
  sigma <- initial$sigma
  missing <- is.na(x)
  pattern_key <- apply(missing, 1L, paste0, collapse = "")
  pattern_rows <- split(seq_len(n), pattern_key)

  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(maxiter)) {
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
            paste0(
              "the observed-data covariance matrix became singular during ",
              "normal-theory EM estimation"
            ),
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
      stop(
        paste0(
          "normal-theory maximum-likelihood estimation produced non-finite ",
          "values"
        ),
        call. = FALSE
      )
    }

    if (difference <= tolerance * parameter_scale) {
      converged <- TRUE
      break
    }
  }

  if (!converged) {
    warning(
      sprintf(
        paste0(
          "the EM algorithm for normal-theory maximum-likelihood estimation ",
          "did not converge within %d iterations"
        ),
        maxiter
      ),
      call. = FALSE
    )
  }

  eigenvalues <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(1, max(abs(eigenvalues)))
  if (any(eigenvalues <= sqrt(.Machine$double.eps) * scale)) {
    stop(
      paste0(
        "the normal-theory maximum-likelihood covariance matrix is singular ",
        "or nearly singular"
      ),
      call. = FALSE
    )
  }

  list(mu = mu, sigma = sigma, iterations = iter, converged = converged)
}


check_mcar_matrix <- function(x, test_label) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop(
        sprintf("all variables supplied to %s must be numeric", test_label),
        call. = FALSE
      )
    }
    x <- as.matrix(x)
  } else if (is.matrix(x)) {
    if (!is.numeric(x)) {
      stop(
        sprintf("the matrix supplied to %s must be numeric", test_label),
        call. = FALSE
      )
    }
  } else {
    stop(
      paste0(
        "'object' must be a fitted 'geer' object, numeric matrix, or numeric ",
        "data frame"
      ),
      call. = FALSE
    )
  }
  storage.mode(x) <- "double"
  if (ncol(x) < 2L) {
    stop(
      sprintf(
        "%s requires at least two variables or repeated measurements",
        test_label
      ),
      call. = FALSE
    )
  }
  if (any(is.infinite(x))) {
    stop(
      sprintf(
        "data supplied to %s must not contain infinite values",
        test_label
      ),
      call. = FALSE
    )
  }
  original_nrow <- nrow(x)
  row_index <- seq_len(original_nrow)
  all_missing <- rowSums(!is.na(x)) == 0L
  if (any(all_missing)) {
    warning(
      sprintf(
        "%d row(s) with no observed values were omitted",
        sum(all_missing)
      ),
      call. = FALSE
    )
    x <- x[!all_missing, , drop = FALSE]
    row_index <- row_index[!all_missing]
  }
  if (nrow(x) < 2L) {
    stop(
      sprintf("%s requires at least two rows with observed values", test_label),
      call. = FALSE
    )
  }
  observed_no <- colSums(!is.na(x))
  labels <- colnames(x)
  if (is.null(labels)) {
    labels <- as.character(seq_len(ncol(x)))
  }
  if (any(observed_no < 2L)) {
    stop(
      sprintf(
        "%s requires at least two observed values for variable(s) %s",
        test_label,
        paste(labels[observed_no < 2L], collapse = ", ")
      ),
      call. = FALSE
    )
  }
  attr(x, "row.index") <- row_index
  attr(x, "original.nrow") <- original_nrow
  x
}


get_mcar_data <- function(object, data) {
  if (is.null(data)) {
    data <- object$data
  }
  if (is.null(data)) {
    stop(
      paste0(
        "the original data are not available in the fitted object; supply ",
        "them through 'data'"
      ),
      call. = FALSE
    )
  }
  data
}


build_mcar_model_frame <- function(mcall, eval_env, formula_env, what) {
  tryCatch(
    eval(mcall, envir = eval_env, enclos = formula_env),
    error = function(e) {
      stop(
        sprintf("could not %s: %s", what, conditionMessage(e)),
        call. = FALSE
      )
    }
  )
}


check_mcar_univariate <- function(response) {
  if (is.matrix(response) || is.data.frame(response)) {
    stop("the MCAR tests require a univariate response", call. = FALSE)
  }
  invisible(response)
}


extract_mcar_matrix <- function(object, data = NULL) {
  data <- get_mcar_data(object, data)

  mcall <- object$call
  keep <- c("formula", "data", "subset", "id", "repeated", "na.action")
  mcall <- mcall[c(1L, match(keep, names(mcall), 0L))]
  mcall[[1L]] <- quote(stats::model.frame)
  mcall$formula <- stats::formula(object$terms)
  mcall$data <- quote(.geer_mcar_data)
  mcall$na.action <- quote(stats::na.pass)
  mcall$drop.unused.levels <- TRUE

  formula_env <- environment(stats::formula(object$terms))
  if (is.null(formula_env)) {
    formula_env <- parent.frame()
  }
  eval_env <- new.env(parent = formula_env)
  assign(".geer_mcar_data", data, envir = eval_env)
  model_frame <- build_mcar_model_frame(
    mcall,
    eval_env,
    formula_env,
    "reconstruct the data used by the fit"
  )

  response <- stats::model.response(model_frame)
  check_mcar_univariate(response)
  if (!is.numeric(response)) {
    stop(
      "the MCAR tests require a numeric response",
      call. = FALSE
    )
  }

  id_repeated <- extract_geer_id_repeated(
    model_frame,
    length(response),
    warn_character = FALSE
  )
  id <- id_repeated$id
  repeated <- id_repeated$repeated
  repeated_raw <- stats::model.extract(model_frame, "repeated")
  repeated_labels <- if (is.null(repeated_raw)) {
    as.character(seq_len(max(repeated)))
  } else {
    levels(factor(repeated_raw))
  }

  response_matrix <- matrix(
    NA_real_,
    nrow = max(id),
    ncol = max(repeated),
    dimnames = list(
      levels(factor(stats::model.extract(model_frame, "id"))),
      repeated_labels
    )
  )
  response_matrix[cbind(id, repeated)] <- as.numeric(response)
  response_matrix
}
