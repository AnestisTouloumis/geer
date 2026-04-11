#' @title emmeans Support for `geer` Objects
#'
#' @description
#' Methods that allow fitted `geer` objects to work with the
#' \pkg{emmeans} package.
#'
#' These functions are used internally by \pkg{emmeans} and are not usually
#' called directly by end users. Once both packages are installed, functions
#' such as `emmeans::emmeans()` and `emmeans::ref_grid()` can be applied to
#' fitted `geer` objects.
#'
#' By default, `emmeans` calculations for `geer` objects use the robust
#' covariance matrix. Alternative covariance estimators may be requested via
#' `vcov.method`, with supported character values
#' `"robust"`, `"bias-corrected"`, `"df-adjusted"`, and `"naive"`.
#' The alias `"model"` is also accepted for `"naive"`.
#'
#' In line with the large-sample inference used elsewhere in \pkg{geer},
#' these methods use asymptotic inference with degrees of freedom equal to
#' `Inf`.
#'
#' @param object A fitted object of class `"geer"`.
#' @param data Optional data frame used by `emmeans` to recover the predictors.
#' @param trms The terms component supplied by `emmeans`.
#' @param xlev The factor levels supplied by `emmeans`.
#' @param grid The reference grid supplied by `emmeans`.
#' @param vcov.method Covariance specification to use for `emmeans`
#'   calculations. This may be a character string naming one of the supported
#'   covariance estimators, or a covariance matrix or function accepted by
#'   `emmeans` through its `vcov.` mechanism.
#' @param cov_type Optional alias for `vcov.method`.
#' @param vcov. Optional covariance matrix or function supplied through the
#'   standard `emmeans` `vcov.` argument. When provided, it takes precedence
#'   over `vcov.method`.
#' @param misc Optional list passed through by `emmeans`.
#' @param ... Additional arguments passed through from \pkg{emmeans}.
#'
#' @return
#' `recover_data.geer()` returns a recovered predictor data set for use by
#' \pkg{emmeans}. `emm_basis.geer()` returns the basis components required by
#' \pkg{emmeans} to construct an `emmGrid` object.
#'
#' @examples
#' if (requireNamespace("emmeans", quietly = TRUE)) {
#'   data("cerebrovascular", package = "geer")
#'
#'   fit <- geewa_binary(
#'     formula = ecg ~ treatment + factor(period),
#'     id = id,
#'     data = cerebrovascular,
#'     link = "logit",
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
