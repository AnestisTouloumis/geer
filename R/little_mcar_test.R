validate_little_mcar_matrix <- function(x) {
  if (is.data.frame(x)) {
    numeric_cols <- vapply(x, is.numeric, logical(1))
    if (!all(numeric_cols)) {
      stop("all variables supplied to Little's MCAR test must be numeric", call. = FALSE)
    }
    x <- as.matrix(x)
  } else if (is.matrix(x)) {
    if (!is.numeric(x)) {
      stop("the matrix supplied to Little's MCAR test must be numeric", call. = FALSE)
    }
  } else {
    stop("'object' must be a fitted 'geer' object, numeric matrix, or numeric data frame", call. = FALSE)
  }

  storage.mode(x) <- "double"
  if (nrow(x) < 2L) {
    stop("Little's MCAR test requires at least two rows", call. = FALSE)
  }
  if (ncol(x) < 2L) {
    stop("Little's MCAR test requires at least two variables or repeated measurements", call. = FALSE)
  }
  if (any(is.infinite(x))) {
    stop("data supplied to Little's MCAR test cannot contain infinite values", call. = FALSE)
  }
  observed_no <- colSums(!is.na(x))
  labels <- colnames(x)
  if (is.null(labels)) labels <- as.character(seq_len(ncol(x)))
  if (any(observed_no < 2L)) {
    bad <- which(observed_no < 2L)
    stop(
      sprintf(
        "Little's MCAR test requires at least two observed values for variable(s) %s",
        paste(labels[bad], collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (ncol(x) > 1L) {
    pairs <- utils::combn(ncol(x), 2L)
    overlap <- vapply(
      seq_len(ncol(pairs)),
      function(j) {
        sum(!is.na(x[, pairs[1L, j]]) & !is.na(x[, pairs[2L, j]]))
      },
      integer(1)
    )
    if (any(overlap < 1L)) {
      bad_pair <- pairs[, which(overlap < 1L)[1L]]
      stop(
        sprintf(
          paste0(
            "Little's MCAR test requires at least one jointly observed value ",
            "for every variable pair; '%s' and '%s' do not meet this requirement"
          ),
          labels[bad_pair[1L]],
          labels[bad_pair[2L]]
        ),
        call. = FALSE
      )
    }
  }
  x
}

little_mcar_initial_parameters <- function(x) {
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

little_mcar_em <- function(x, maxit, tol) {
  n <- nrow(x)
  p <- ncol(x)
  initial <- little_mcar_initial_parameters(x)
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

little_mcar_bivariate_exact <- function(x) {
  n <- nrow(x)
  if (ncol(x) != 2L || n <= 2L) return(NULL)

  observed_no <- colSums(!is.na(x))
  complete_col <- which(observed_no == n)
  incomplete_col <- which(observed_no < n)
  if (length(complete_col) != 1L || length(incomplete_col) != 1L) return(NULL)

  missing_group <- is.na(x[, incomplete_col])
  if (!any(missing_group) || all(missing_group)) return(NULL)

  y <- x[, complete_col]
  group_mean_observed <- mean(y[!missing_group])
  group_mean_missing <- mean(y[missing_group])
  overall_mean <- mean(y)
  n_observed <- sum(!missing_group)
  n_missing <- sum(missing_group)

  ss_between <- n_observed * (group_mean_observed - overall_mean)^2 +
    n_missing * (group_mean_missing - overall_mean)^2
  ss_within <- sum((y[!missing_group] - group_mean_observed)^2) +
    sum((y[missing_group] - group_mean_missing)^2)
  df2 <- n - 2L
  f_statistic <- if (ss_within == 0) Inf else ss_between / (ss_within / df2)
  p_value <- stats::pf(f_statistic, df1 = 1, df2 = df2, lower.tail = FALSE)

  # Little (1988), Section 3.4: d^2 = (n - 1) F / (n - 2 + F).
  implied_statistic <- if (is.infinite(f_statistic)) {
    n - 1
  } else {
    (n - 1) * f_statistic / (n - 2 + f_statistic)
  }

  list(
    statistic = f_statistic,
    df1 = 1L,
    df2 = as.integer(df2),
    p_value = p_value,
    implied_little_statistic = implied_statistic
  )
}

little_mcar_calculate <- function(x, maxit, tol) {
  missing <- is.na(x)
  n <- nrow(x)
  p <- ncol(x)
  pattern_key <- apply(missing, 1L, paste0, collapse = "")
  pattern_rows <- split(seq_len(n), pattern_key)
  missing_patterns <- length(pattern_rows)

  if (!any(missing)) {
    sigma_ml <- stats::cov(x) * (n - 1) / n
    sigma_corrected <- stats::cov(x)
    return(list(
      statistic = 0,
      df = 0L,
      p_value = 1,
      missing_patterns = 1L,
      iterations = 0L,
      converged = TRUE,
      mu = colMeans(x),
      sigma = sigma_corrected,
      sigma_ml = sigma_ml,
      bivariate_exact = NULL
    ))
  }

  em <- little_mcar_em(x, maxit = maxit, tol = tol)
  sigma_corrected <- (n / (n - 1)) * em$sigma
  statistic <- 0
  df_terms <- 0L

  for (rows in pattern_rows) {
    observed <- which(!missing[rows[1L], ])
    k <- length(observed)
    df_terms <- df_terms + k
    if (!k) next

    pattern_mean <- colMeans(x[rows, observed, drop = FALSE])
    difference <- pattern_mean - em$mu[observed]
    sigma_observed <- sigma_corrected[observed, observed, drop = FALSE]
    solved <- tryCatch(
      solve(sigma_observed, difference),
      error = function(e) {
        stop(
          "a covariance submatrix required for Little's MCAR statistic is singular",
          call. = FALSE
        )
      }
    )
    statistic <- statistic + length(rows) * as.numeric(crossprod(difference, solved))
  }

  df <- as.integer(df_terms - p)
  if (df <= 0L) {
    stop(
      "Little's MCAR test has non-positive degrees of freedom for this missing-data pattern",
      call. = FALSE
    )
  }

  statistic <- max(as.numeric(statistic), 0)
  asymptotic_p_value <- stats::pchisq(statistic, df = df, lower.tail = FALSE)
  bivariate_exact <- little_mcar_bivariate_exact(x)

  if (!is.null(bivariate_exact) &&
      !isTRUE(all.equal(
        statistic,
        bivariate_exact$implied_little_statistic,
        tolerance = 1e-6,
        check.attributes = FALSE
      ))) {
    warning(
      paste0(
        "the bivariate monotone small-sample identity in Little (1988) was ",
        "not reproduced to numerical tolerance; using the asymptotic reference"
      ),
      call. = FALSE
    )
    bivariate_exact <- NULL
  }

  list(
    statistic = statistic,
    df = df,
    p_value = asymptotic_p_value,
    missing_patterns = as.integer(missing_patterns),
    iterations = as.integer(em$iterations),
    converged = isTRUE(em$converged),
    mu = em$mu,
    sigma = sigma_corrected,
    sigma_ml = em$sigma,
    bivariate_exact = bivariate_exact
  )
}

little_mcar_has_binary_variable <- function(x) {
  any(vapply(
    seq_len(ncol(x)),
    function(j) {
      values <- unique(x[!is.na(x[, j]), j])
      length(values) >= 2L && all(values %in% c(0, 1))
    },
    logical(1)
  ))
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
  mcall$formula <- formula(object$terms)
  mcall$data <- quote(.geer_mcar_data)
  mcall$na.action <- quote(stats::na.pass)
  mcall$drop.unused.levels <- TRUE

  formula_env <- environment(formula(object$terms))
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

  response <- model.response(model_frame)
  if (is.matrix(response) || is.data.frame(response)) {
    stop(
      "Little's MCAR test for a 'geer' object currently requires a univariate response at each repeated measurement",
      call. = FALSE
    )
  }
  if (!is.numeric(response)) {
    stop("the response used for Little's MCAR test must be numeric", call. = FALSE)
  }

  id_raw <- model.extract(model_frame, "id")
  if (is.null(id_raw)) stop("'id' could not be recovered from the fitted model", call. = FALSE)
  if (anyNA(id_raw)) stop("'id' cannot contain missing values for Little's MCAR test", call. = FALSE)
  id <- as.numeric(factor(id_raw))

  repeated_raw <- model.extract(model_frame, "repeated")
  if (is.null(repeated_raw)) {
    repeated <- ave(id, id, FUN = seq_along)
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


#' Little's Test for Missing Completely at Random Data
#'
#' @description
#' Performs Little's (1988) test of the null hypothesis that
#' missing values are missing completely at random (MCAR). For a fitted
#' \code{geer} object, the repeated response is reconstructed in wide form
#' from the original data before rows with missing responses were omitted by
#' model fitting.
#'
#' @param object a fitted object of class \code{"geer"}, a numeric matrix, or
#'   a numeric data frame. For matrix or data-frame input, rows are independent
#'   units and columns are variables or repeated measurements.
#' @param data optional original data used to fit \code{object}. This is only
#'   used when \code{object} is a \code{geer} fit and is useful when the
#'   original data cannot be recovered from the fitted object.
#' @param maxit positive integer giving the maximum number of EM iterations
#'   used to obtain the multivariate-normal maximum-likelihood estimates.
#'   Defaults to 1000.
#' @param tol positive convergence tolerance for the EM algorithm. Defaults to
#'   \code{1e-08}.
#' @param reference reference distribution used for the p-value. \code{"auto"}
#'   uses Little's exact bivariate monotone result when it applies and otherwise
#'   uses the asymptotic chi-squared distribution. \code{"asymptotic"}
#'   always uses the chi-squared reference.
#'
#' @details
#' Let the rows be divided into missing-data patterns. For pattern \eqn{r},
#' let \eqn{n_r} be its number of rows and \eqn{\bar y_r} the mean vector
#' of the variables observed in that pattern. Let \eqn{\hat\mu} and
#' \eqn{\hat\Sigma} be the common mean vector and covariance matrix estimated
#' by maximum likelihood under an ignorable multivariate-normal missing-data
#' model, and define Little's degrees-of-freedom correction
#' \eqn{\widetilde\Sigma = n\hat\Sigma/(n-1)}. Little's statistic is
#' \deqn{d^2 = \sum_r n_r (\bar y_r-\hat\mu_r)^T
#' \widetilde\Sigma_r^{-1}(\bar y_r-\hat\mu_r).}
#' Under MCAR, the statistic is asymptotically chi-squared with
#' \deqn{\mathrm{df} = \sum_r p_r - p,}
#' where \eqn{p_r} is the number of variables observed in pattern \eqn{r}
#' and \eqn{p} is the total number of variables.
#'
#' The common mean and covariance are estimated internally with an EM
#' algorithm, so no additional package is required. Following Little (1988),
#' the multivariate-normal maximum-likelihood covariance estimate
#' \eqn{\hat\Sigma} is multiplied by \eqn{n/(n-1)} before its pattern-specific
#' submatrices are used in the test statistic.
#'
#' For two variables with one variable observed for every row and missingness
#' confined to the other variable, Little shows that the small-sample null
#' distribution can be obtained from an ordinary two-group ANOVA (equivalently,
#' a pooled two-sample t test):
#' \deqn{d^2 = \frac{(n-1)F}{n-2+F}, \qquad F \sim F_{1,n-2}.}
#' With \code{reference = "auto"}, this exact normal-theory reference is
#' used automatically in that special bivariate monotone case. For all other
#' missing-data patterns, the large-sample chi-squared reference is used. Little
#' also derives a more general small-sample distribution for monotone patterns
#' as a sum of transformed independent F variables; that nonstandard reference
#' distribution is not evaluated here.
#'
#' For a \code{geer} fit, only the repeated response is tested. The function
#' reconstructs a subject-by-repeated-measure matrix from \code{id} and
#' \code{repeated}. If \code{repeated} was omitted during fitting, the
#' within-cluster row order in the original data defines the repeated-measure
#' positions. Consequently, an entirely absent row cannot be distinguished
#' from a measurement occasion that never existed unless \code{repeated} was
#' supplied explicitly.
#'
#' Little's test is based on a multivariate-normal working model for the
#' variables being tested. Little (1988) notes that the test is most appropriate
#' for quantitative variables and recommends contingency-table methods for
#' categorical variables. The function therefore warns when a binary variable
#' is detected. The paper also reports that the asymptotic chi-squared test can
#' be conservative in small samples.
#'
#' As emphasized by Hardin and Hilbe (2013), this is a diagnostic for the
#' missingness mechanism rather than a test of the fitted GEE mean model itself.
#' A small p-value provides evidence against MCAR; a large p-value does not
#' establish that MCAR holds.
#'
#' @return
#' An object of class \code{"htest"}. In addition to the standard components,
#' the object contains:
#' \item{missing.patterns}{the number of distinct missing-data patterns.}
#' \item{n}{the number of independent rows or clusters tested.}
#' \item{variables}{the number of variables or repeated measurements tested.}
#' \item{iterations}{the number of EM iterations used.}
#' \item{converged}{whether the EM algorithm converged.}
#' \item{mean}{the maximum-likelihood estimate of the common mean vector.}
#' \item{covariance}{the degrees-of-freedom-corrected covariance matrix used
#' in Little's statistic.}
#' \item{ml.covariance}{the uncorrected multivariate-normal maximum-likelihood
#' covariance estimate.}
#' \item{reference}{the reference distribution used for the reported p-value.}
#' \item{asymptotic.df}{the degrees of freedom for the large-sample
#' chi-squared reference.}
#' \item{asymptotic.p.value}{the large-sample chi-squared p-value.}
#' \item{exact.p.value}{the exact bivariate monotone p-value when available,
#' otherwise \code{NULL}.}
#' \item{exact.f.statistic}{the corresponding bivariate ANOVA F statistic
#' when available, otherwise \code{NULL}.}
#'
#' @references
#' Little, R.J.A. (1988). A test of missing completely at random for
#' multivariate data with missing values. \emph{Journal of the American
#' Statistical Association}, \bold{83}, 1198--1202.
#' \doi{10.1080/01621459.1988.10478722}
#'
#' Hardin, J.W. and Hilbe, J.M. (2013). \emph{Generalized Estimating
#' Equations}, 2nd ed. Chapman and Hall/CRC, Section 4.6.
#'
#' @seealso \code{\link{runs_test}}, \code{\link{residuals.geer}},
#'   \code{\link{geewa}}, \code{\link{geewa_binary}}.
#'
#' @examples
#' x <- data.frame(
#'   visit1 = c(2.1, 4.3, 3.2, 5.8, 7.1, 6.4, 8.2, 9.5, 10.1, 11.3),
#'   visit2 = c(1.5, 3.8, 2.9, 6.2, NA, 5.1, 8.9, NA, 9.4, 12.0),
#'   visit3 = c(3.0, 4.9, NA, 7.1, 6.3, 7.8, NA, 10.2, 11.5, 13.1)
#' )
#' little_mcar_test(x)
#'
#' @export
little_mcar_test <- function(object, data = NULL, maxit = 1000L, tol = 1e-8,
                             reference = c("auto", "asymptotic")) {
  if (!is.numeric(maxit) || length(maxit) != 1L || is.na(maxit) ||
      !is.finite(maxit) || maxit < 1 || maxit > .Machine$integer.max ||
      maxit != floor(maxit)) {
    stop("'maxit' must be a positive integer", call. = FALSE)
  }
  maxit <- as.integer(maxit)
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
      !is.finite(tol) || tol <= 0) {
    stop("'tol' must be a single positive finite number", call. = FALSE)
  }
  reference <- match.arg(reference)

  data_name <- deparse(substitute(object))
  if (inherits(object, "geer")) {
    object <- check_geer_object(object)
    x <- extract_geer_mcar_matrix(object, data = data)
    data_name <- "repeated response from fitted geer object"
  } else {
    if (!is.null(data)) {
      stop("'data' is only used when 'object' is a fitted 'geer' object", call. = FALSE)
    }
    x <- object
  }

  x <- validate_little_mcar_matrix(x)
  if (anyNA(x) && little_mcar_has_binary_variable(x)) {
    warning(
      paste0(
        "Little (1988) notes that this MCAR test is most appropriate for ",
        "quantitative variables; consider contingency-table methods for ",
        "categorical variables"
      ),
      call. = FALSE
    )
  }

  result <- little_mcar_calculate(x, maxit = maxit, tol = tol)
  exact <- result$bivariate_exact
  use_exact <- identical(reference, "auto") && !is.null(exact)

  if (use_exact) {
    statistic <- c("d^2" = result$statistic)
    parameter <- c(df1 = exact$df1, df2 = exact$df2)
    p_value <- exact$p_value
    method <- "Little's MCAR test (exact bivariate monotone reference)"
    reference_used <- "exact bivariate F"
  } else {
    statistic <- c("X-squared" = result$statistic)
    parameter <- c(df = result$df)
    p_value <- result$p_value
    method <- "Little's MCAR test"
    reference_used <- "asymptotic chi-squared"
  }

  structure(
    list(
      statistic = statistic,
      parameter = parameter,
      p.value = p_value,
      method = method,
      data.name = data_name,
      alternative = "data are not missing completely at random",
      missing.patterns = result$missing_patterns,
      n = nrow(x),
      variables = ncol(x),
      iterations = result$iterations,
      converged = result$converged,
      mean = result$mu,
      covariance = result$sigma,
      ml.covariance = result$sigma_ml,
      reference = reference_used,
      asymptotic.df = result$df,
      asymptotic.p.value = result$p_value,
      exact.p.value = if (is.null(exact)) NULL else exact$p_value,
      exact.f.statistic = if (is.null(exact)) NULL else exact$statistic
    ),
    class = "htest"
  )
}
