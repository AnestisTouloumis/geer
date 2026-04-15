#' @title
#' emmeans Support for geer Objects
#'
#' @description
#' Methods that allow fitted \code{geer} objects to work with the
#' \pkg{emmeans} package.
#'
#' These functions are used internally by \pkg{emmeans} and are not usually
#' called directly by end users. Once both packages are installed, functions
#' such as \code{emmeans::emmeans()} and \code{emmeans::ref_grid()} can be
#' applied to fitted \code{geer} objects.
#'
#' By default, \pkg{emmeans} calculations for \code{geer} objects use the
#' robust covariance matrix. Alternative covariance estimators may be
#' requested via \code{vcov.method}, with supported character values
#' \code{"robust"}, \code{"bias-corrected"}, \code{"df-adjusted"}, and
#' \code{"naive"}. The alias \code{"model"} is also accepted for
#' \code{"naive"}.
#'
#' In line with the large-sample inference used elsewhere in \pkg{geer},
#' these methods use asymptotic inference with degrees of freedom equal to
#' \code{Inf}.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param data optional data frame used by \pkg{emmeans} to recover the
#'        model predictors.
#' @param trms the terms component supplied by \pkg{emmeans}.
#' @param xlev the factor levels supplied by \pkg{emmeans}.
#' @param grid the reference grid supplied by \pkg{emmeans}.
#' @param vcov.method covariance specification to use for \pkg{emmeans}
#'        calculations. This may be a character string specifying one of the
#'        supported covariance estimators (\code{"robust"},
#'        \code{"bias-corrected"}, \code{"df-adjusted"}, \code{"naive"}, or
#'        the alias \code{"model"} for \code{"naive"}), or a covariance matrix
#'        or function accepted by \pkg{emmeans} through its \code{vcov.}
#'        mechanism. Defaults to \code{"robust"}.
#' @param cov_type optional alias for \code{vcov.method}.
#' @param vcov. optional covariance matrix or function supplied through the
#'        standard \pkg{emmeans} \code{vcov.} argument. When provided, it
#'        takes precedence over \code{vcov.method}.
#' @param misc optional list passed through by \pkg{emmeans}.
#' @param ... additional arguments passed through from \pkg{emmeans}.
#'
#' @return
#' \code{recover_data.geer()} returns a recovered predictor data frame for use
#' by \pkg{emmeans}. \code{emm_basis.geer()} returns the basis components
#' required by \pkg{emmeans} to construct an \code{emmGrid} object.
#'
#' @seealso \code{\link{vcov.geer}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
#'
#' @examples
#' if (requireNamespace("emmeans", quietly = TRUE)) {
#'   data("cerebrovascular", package = "geer")
#'
#'   fit <- geewa_binary(
#'     formula = ecg ~ treatment + factor(period),
#'     link = "logit",
#'     data = cerebrovascular,
#'     id = id,
#'     orstr = "exchangeable"
#'   )
#'
#'   emmeans::emmeans(fit, ~ treatment)
#'   emmeans::emmeans(fit, ~ treatment, type = "response")
#'   emmeans::emmeans(fit, ~ treatment, vcov.method = "naive")
#' }
#'
#' @name emmeans-support
NULL

normalize_emmeans_vcov_method <- function(vcov.method) {
  if (is.null(vcov.method)) {
    return("robust")
  }
  if (!is.character(vcov.method) || length(vcov.method) != 1L || is.na(vcov.method)) {
    stop(
      "'vcov.method' must be a single character value, matrix, or function",
      call. = FALSE
    )
  }
  key <- tolower(vcov.method)
  key <- gsub("[._]", "-", key)
  choices <- c("robust", "bias-corrected", "df-adjusted", "naive", "model")
  idx <- pmatch(key, choices)
  if (is.na(idx)) {
    stop(
      "'vcov.method' must identify one of: robust, bias-corrected, df-adjusted, naive, model",
      call. = FALSE
    )
  }
  choice <- choices[[idx]]
  if (identical(choice, "model")) {
    choice <- "naive"
  }
  choice
}


resolve_emmeans_vcov <- function(object,
                                 vcov.method = "robust",
                                 vcov. = NULL,
                                 ...) {
  if (!is.null(vcov.)) {
    return(emmeans::.my.vcov(object, vcov. = vcov., ...))
  }
  if (is.character(vcov.method) || is.null(vcov.method)) {
    cov_type <- normalize_emmeans_vcov_method(vcov.method)
    return(vcov(object, cov_type = cov_type))
  }
  if (is.matrix(vcov.method) || is.function(vcov.method)) {
    return(emmeans::.my.vcov(object, vcov. = vcov.method, ...))
  }
  stop(
    "'vcov.method' must be a character value, covariance matrix, or function",
    call. = FALSE
  )
}


#' @rdname emmeans-support
#' @export
recover_data.geer <- function(object,
                              data = NULL,
                              ...) {
  object <- check_geer_object(object)

  recovered_data <- if (is.null(data) && is.data.frame(object$data)) {
    object$data
  } else {
    data
  }
  emmeans::recover_data(
    object$call,
    trms = delete.response(object$terms),
    na.action = object$na.action,
    data = recovered_data,
    pwts = object$prior.weights,
    ...
  )
}


#' @rdname emmeans-support
#' @export
emm_basis.geer <- function(object,
                           trms,
                           xlev,
                           grid,
                           vcov.method = "robust",
                           cov_type = NULL,
                           vcov. = NULL,
                           misc = NULL,
                           ...) {
  object <- check_geer_object(object)

  if (!is.null(cov_type)) {
    vcov.method <- cov_type
  }
  model_frame <- model.frame(trms, grid, na.action = na.pass, xlev = xlev)
  design_matrix <- model.matrix(trms, model_frame, contrasts.arg = object$contrasts)
  coefficient_vector <- coef(object)
  coefficient_names <- names(coefficient_vector)
  if (is.null(coefficient_names) || !length(coefficient_names)) {
    stop("'coefficients' must be a named numeric vector", call. = FALSE)
  }
  if (is.null(colnames(design_matrix))) {
    stop("the model matrix reconstructed for 'emmeans' must have column names", call. = FALSE)
  }
  missing_cols <- setdiff(coefficient_names, colnames(design_matrix))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "failed to reconstruct the model matrix for 'emmeans'; missing columns: %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  design_matrix <- design_matrix[, coefficient_names, drop = FALSE]
  vcov_matrix <- resolve_emmeans_vcov(
    object,
    vcov.method = vcov.method,
    vcov. = vcov.,
    ...
  )
  vcov_matrix <- vcov_matrix[coefficient_names, coefficient_names, drop = FALSE]
  nbasis <- if (anyNA(coefficient_vector)) {
    estimability::nonest.basis(object$qr)
  } else {
    estimability::all.estble
  }
  dfargs <- list(df = Inf)
  dffun <- function(k, dfargs) dfargs$df
  misc <- misc %||% list()
  misc <- emmeans::.std.link.labels(object$family, misc)
  if (identical(object$family$family, "gaussian") &&
      is.numeric(object$phi) && length(object$phi) == 1L && is.finite(object$phi)) {
    misc$sigma <- sqrt(object$phi)
  }
  offset <- emmeans::.get.offset(trms, grid)
  ans <- list(
    X = design_matrix,
    bhat = as.numeric(coefficient_vector),
    nbasis = nbasis,
    V = vcov_matrix,
    dffun = dffun,
    dfargs = dfargs,
    misc = misc
  )
  if (!is.null(offset)) {
    ans$offset <- offset
  }
  ans
}
