#' @title
#' emmeans Support for geer Objects
#'
#' @description
#' Methods that allow \pkg{emmeans} functions, such as
#' \code{emmeans::emmeans()}, to be applied to fitted \code{geer} objects.
#' They are called by \pkg{emmeans} and not usually by the user. Inference is
#' asymptotic, with infinite degrees of freedom, and uses the bias-corrected
#' covariance matrix unless \code{vcov.method} requests another estimator.
#' Gaussian models supply their working residual standard deviation through
#' \code{\link{sigma.geer}}.
#' Stored-data recovery uses the rows retained by the fit, including its
#' subsetting and missing-value exclusions.
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param data optional data frame from which \pkg{emmeans} recovers the
#'   predictors. Defaults to \code{NULL}, which uses the data used to fit the
#'   model.
#' @param trms the terms component supplied by \pkg{emmeans}.
#' @param xlev the factor levels supplied by \pkg{emmeans}.
#' @param grid the reference grid supplied by \pkg{emmeans}.
#' @param vcov.method character string specifying the covariance matrix
#'   estimator, one of the estimators described in \code{\link{vcov.geer}},
#'   or a covariance matrix or function accepted by \pkg{emmeans}. Defaults to
#'   \code{"bias-corrected"}.
#' @param cov_type optional alias of \code{vcov.method}, which it overrides.
#'   Defaults to \code{NULL}.
#' @param vcov. optional covariance matrix or function supplied through the
#'   \pkg{emmeans} argument \code{vcov.}, which overrides \code{cov_type} and
#'   \code{vcov.method}. Defaults to \code{NULL}.
#' @param misc optional list passed through by \pkg{emmeans}. Defaults to
#'   \code{NULL}.
#' @param ... additional framework arguments, forwarded to the recovery or
#'   covariance machinery where applicable. Arguments unused by those routes
#'   are ignored for compatibility with \pkg{emmeans}.
#'
#' @return
#' \code{recover_data.geer()} returns the predictor data, and
#' \code{emm_basis.geer()} the components from which \pkg{emmeans} builds a
#' reference grid.
#'
#' @seealso \code{\link{vcov.geer}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
#'
#' @examples
#' if (requireNamespace("emmeans", quietly = TRUE)) {
#'   data("cerebrovascular", package = "geer")
#'   fit <- geewa_binary(
#'     formula = ecg ~ treatment + factor(period),
#'     link = "logit",
#'     data = cerebrovascular,
#'     id = id,
#'     repeated = period,
#'     orstr = "exchangeable"
#'   )
#'   print(emmeans::emmeans(fit, specs = ~ treatment))
#'   print(emmeans::emmeans(fit, specs = ~ treatment, type = "response"))
#' }
#'
#' @name emmeans-support
NULL


normalize_emmeans_vcov_method <- function(vcov.method) {
  if (is.null(vcov.method)) {
    return("bias-corrected")
  }
  if (!is.character(vcov.method) || length(vcov.method) != 1L ||
        is.na(vcov.method)) {
    stop(
      paste0(
        "'vcov.method' must be a character string, a covariance matrix, or a ",
        "function"
      ),
      call. = FALSE
    )
  }
  match_choice_with_default(
    x = vcov.method,
    choices = geer_cov_type_choices,
    name = "vcov.method"
  )
}


get_emmeans_vcov <- function(object,
                             vcov.method = "bias-corrected",
                             vcov. = NULL,
                             ...) {
  if (!is.null(vcov.)) {
    return(emmeans::.my.vcov(object, vcov. = vcov., ...))
  }
  if (is.character(vcov.method) || is.null(vcov.method)) {
    cov_type <- normalize_emmeans_vcov_method(vcov.method)
    return(stats::vcov(object, cov_type = cov_type))
  }
  if (is.matrix(vcov.method) || is.function(vcov.method)) {
    return(emmeans::.my.vcov(object, vcov. = vcov.method, ...))
  }
  stop(
    paste0(
      "'vcov.method' must be a character string, a covariance matrix, or a ",
      "function"
    ),
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
    rows <- object$frame_rows
    if (is.null(rows)) {
      model_call <- object$call
      model_call$data <- object$data
      model_call$formula <- object$formula
      frame <- build_geer_model_frame(
        model_call, env = environment(object$formula)
      )
      rows <- get_geer_frame_rows(frame, object$data)
    }
    if (!is.numeric(rows) || length(rows) != nrow(object$x) ||
          any(!is.finite(rows)) || any(rows != floor(rows)) ||
          any(rows < 1L | rows > nrow(object$data)) || anyDuplicated(rows)) {
      stop(
        "recover_data.geer failed: the fitting rows could not be recovered",
        call. = FALSE
      )
    }
    object$data[rows, , drop = FALSE]
  } else {
    data
  }
  prior_weights <- object$prior.weights
  if (!is.null(object$row_order) &&
        length(object$row_order) == length(prior_weights)) {
    prior_weights <- prior_weights[order(object$row_order)]
  }
  emmeans::recover_data(
    object$call,
    trms = stats::delete.response(object$terms),
    na.action = object$na.action,
    data = recovered_data,
    pwts = prior_weights,
    ...
  )
}


#' @rdname emmeans-support
#' @export
emm_basis.geer <- function(object,
                           trms,
                           xlev,
                           grid,
                           vcov.method = "bias-corrected",
                           cov_type = NULL,
                           vcov. = NULL,
                           misc = NULL,
                           ...) {
  object <- check_geer_object(object)
  if (!is.null(cov_type)) {
    vcov.method <- match_cov_type(cov_type)
  }
  model_frame <- stats::model.frame(
    trms,
    grid,
    na.action = stats::na.pass,
    xlev = xlev
  )
  design_matrix <- stats::model.matrix(
    trms,
    model_frame,
    contrasts.arg = object$contrasts
  )
  coefficient_vector <- stats::coef(object)
  coefficient_names <- names(coefficient_vector)
  if (is.null(coefficient_names) || !length(coefficient_names)) {
    stop("'coefficients' must be a named numeric vector", call. = FALSE)
  }
  if (is.null(colnames(design_matrix))) {
    stop(
      paste0(
        "emm_basis.geer failed: the reconstructed model matrix has no column ",
        "names"
      ),
      call. = FALSE
    )
  }
  missing_cols <- setdiff(coefficient_names, colnames(design_matrix))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        paste0(
          "emm_basis.geer failed: the reconstructed model matrix lacks the ",
          "columns %s"
        ),
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  design_matrix <- design_matrix[, coefficient_names, drop = FALSE]
  vcov_matrix <- get_emmeans_vcov(
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
