#' @title
#' Solving (Adjusted) GEE for Binary Responses
#'
#' @inherit geewa description
#'
#' @inheritParams geewa
#' @param link character specifying the link function. Options include
#'        \code{logit}, \code{probit}, \code{cauchit}, \code{cloglog},
#'        \code{identity}, \code{log}, \code{sqrt}, \code{1/mu^2} and
#'        \code{inverse}. By default, \code{link = "logit"}.
#' @param orstr character specifying the odds ratios structure.
#'        Options include \code{"independence"}, \code{"exchangeable"},
#'        \code{"unstructured"} and \code{"fixed"}. By default,
#'        \code{orstr = "independence"}.
#' @param alpha_vector numerical vector specifying odds ratios structure
#'        when \code{orstr == "fixed"}. It is ignored for all other
#'        possible values of \code{orstr}.
#'
#' @inherit geewa details
#'
#' @inherit geewa return
#'
#' @author Anestis Touloumis
#'
#' @examples
#' data("respiratory")
#' respiratory2 <- respiratory[respiratory$center == "C2", ]
#' fitted_model <- geewa_binary(
#'   formula = status ~ baseline + treatment * gender + visit * age,
#'   id = id, repeated = visit, link = "probit",
#'   data = respiratory2, orstr = "independence",
#'   method = "pgee-jeffreys"
#' )
#' summary(fitted_model, cov_type = "bias-corrected")
#'
#' data("cholecystectomy")
#' fitted_model_gee <- geewa_binary(
#'   formula = pain ~ treatment + gender + age,
#'   id = id, data = cholecystectomy, link = "logit",
#'   orstr = "independence", method = "gee"
#' )
#' summary(fitted_model_gee, cov_type = "bias-corrected")
#' fitted_model_brgee_robust <- update(fitted_model_gee, method = "brgee-robust")
#' summary(fitted_model_brgee_robust, cov_type = "bias-corrected")
#' fitted_model_brgee_naive <- update(fitted_model_gee, method = "brgee-naive")
#' summary(fitted_model_brgee_naive, cov_type = "bias-corrected")
#' fitted_model_brgee_empirical <- update(fitted_model_gee, method = "brgee-empirical")
#' summary(fitted_model_brgee_empirical, cov_type = "bias-corrected")
#'
#' @export
geewa_binary <- function(formula,
                         link = "logit",
                         data = parent.frame(),
                         id,
                         repeated = NULL,
                         control = geer_control(...),
                         orstr = "independence",
                         method = "gee",
                         weights,
                         beta_start = NULL,
                         offset,
                         control_glm = list(...),
                         alpha_vector = NULL,
                         ...) {
  ## call and model frame
  call <- match.call(expand.dots = TRUE)
  mcall <- match.call(expand.dots = FALSE)
  model_frame <- build_geer_model_frame(mcall, env = parent.frame())
  ## link function
  links <- c("logit", "probit", "cauchit", "cloglog", "identity", "log",
             "sqrt", "1/mu^2", "inverse")
  link <- as.character(link)
  if (!isTRUE(link %in% links)) {
    stop("'link' should be one of: logit, probit, cauchit, cloglog, identity, log, sqrt, 1/mu^2, inverse")
  }
  family <- binomial(link = link)
  ## response
  y <- model.response(model_frame, "any")
  if (is.null(y)) stop("response variable not found")
  if (is.factor(y)) {
    y <- as.numeric(y != levels(y)[1L])
  } else if (is.character(y)) {
    yfac <- factor(y)
    y <- as.numeric(yfac != levels(yfac)[1L])
  }
  ## weights
  weights <- as.vector(model.weights(model_frame))
  if (is.null(weights)) {
    weights <- rep.int(1, length(y))
  } else {
    if (!is.numeric(weights)) stop("'weights' should be a numeric vector")
    if (any(!is.finite(weights))) stop("'weights' must be finite")
    if (any(weights <= 0)) stop("'weights' must be strictly positive")
  }
  if (length(weights) != length(y)) stop("'weights' and the response must have the same length")
  ## binomial matrix response: combine trials with user weights
  if (is.matrix(y) && ncol(y) == 2L) {
    trials <- rowSums(y)
    if (any(trials <= 0)) stop("for binomial matrix responses, row sums (trials) must be positive")
    weights <- weights * trials
    y <- y[, 1L] / trials
  }
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop("response variable contains non-finite values")
  ## id (recode to 1..N)
  id_raw <- model.extract(model_frame, "id")
  if (is.null(id_raw)) stop("'id' not found")
  if (anyNA(id_raw)) stop("'id' cannot contain missing values")
  id <- as.numeric(factor(id_raw))
  if (length(id) != length(y)) stop("response variable and 'id' are not of same length")
  ## repeated (recode to 1..T within-cluster if missing)
  repeated <- model.extract(model_frame, "repeated")
  if (is.null(repeated)) {
    repeated <- ave(id, id, FUN = seq_along)
  } else {
    if (anyNA(repeated)) stop("'repeated' cannot contain missing values")
    repeated <- as.numeric(factor(repeated))
  }
  if (length(repeated) != length(y)) stop("response variable and 'repeated' are not of same length")
  if (any(unlist(lapply(split(repeated, id), duplicated)))) {
    stop("'repeated' does not have unique values per 'id'")
  }
  ## offset
  offset <- extract_geer_offset(model_frame, y_length = length(y))
  ## model matrix
  design <- build_geer_design_matrix(model_frame)
  model_terms <- design$terms
  model_matrix <- design$x
  xnames <- design$xnames
  qr_model_matrix <- design$qr
  x_assign <- design$assign
  x_contrasts <- design$contrasts
  ## sort by id then repeated
  ord <- order(id, repeated)
  y <- y[ord]
  model_matrix <- model_matrix[ord, , drop = FALSE]
  weights <- weights[ord]
  offset <- offset[ord]
  id <- id[ord]
  repeated <- repeated[ord]
  attr(model_matrix, "assign") <- x_assign
  attr(model_matrix, "contrasts") <- x_contrasts
  ## control
  control <- normalize_geer_control(control)
  maxiter <- control$maxiter
  tolerance <- control$tolerance
  ## method
  methods <- c("gee",
               "brgee-naive", "brgee-robust", "brgee-empirical",
               "bcgee-naive", "bcgee-robust", "bcgee-empirical",
               "pgee-jeffreys")
  method <- as.character(method)
  if (!isTRUE(method %in% methods)) {
    stop("'method' should be one of: gee, brgee-naive, brgee-robust, brgee-empirical, bcgee-naive, bcgee-robust, bcgee-empirical, pgee-jeffreys")
  }
  ## initial beta
  beta_zero <- compute_start_values_geewa_binary(
    model_matrix = model_matrix,
    y = y,
    weights = weights,
    offset = offset,
    family = family,
    beta_start = beta_start,
    control_glm = control_glm,
    tolerance = tolerance
  )
  ## odds ratios structure
  orstrs <- c("independence", "exchangeable", "unstructured", "fixed")
  if (!isTRUE(orstr %in% orstrs)) {
    stop("'orstr' should be one of: independence, exchangeable, unstructured, fixed")
  }
  if (identical(orstr, "independence")) {
    adding_constant <- NULL
    alpha_vector <- rep.int(1, choose(max(repeated), 2))
  } else if (identical(orstr, "fixed")) {
    adding_constant <- NULL
    pairs_no <- choose(max(repeated), 2)
    if (!is.numeric(alpha_vector)) stop("'alpha_vector' must be a numeric vector")
    if (length(alpha_vector) != pairs_no) stop("'alpha_vector' must be of length ", pairs_no)
    if (any(!is.finite(alpha_vector)) || any(alpha_vector <= 0)) {
      stop("'alpha_vector' must be finite and strictly positive when 'orstr = fixed'")
    }
  } else {
    adding_constant <- control$or_adding
    alpha_vector <- get_marginalized_odds_ratios(
      round(y), id, repeated, weights, adding_constant, orstr
    )
  }
  ## bias-corrected wrapper: solve GEE first
  method_original <- method
  if (method_original %in% c("bcgee-naive", "bcgee-robust", "bcgee-empirical")) {
    method <- "gee"
  }
  ## fit
  geesolver_fit <- fit_bingee_or(
    y, model_matrix, id, repeated, weights, link,
    beta_zero, offset, maxiter, tolerance,
    control$step_maxit, control$step_multi,
    control$jeffreys_power, method, alpha_vector
  )
  ## second pass for bias-corrected estimators
  if (method_original %in% c("bcgee-naive", "bcgee-robust", "bcgee-empirical")) {
    if (geesolver_fit$criterion[ncol(geesolver_fit$beta_mat) - 1L] <= tolerance) {
      method <- sub("bcgee", "brgee", method_original)
      geesolver_fit <- fit_bingee_or(
        y, model_matrix, id, repeated, weights,
        link, as.numeric(geesolver_fit$beta_hat),
        offset, 1, tolerance, 1, 1,
        control$jeffreys_power, method, alpha_vector
      )
      method <- method_original
    } else {
      stop("bias-corrected estimator is undefined because the corresponding GEE model did not converge")
    }
  }
  ## output
  fit <- finalize_geer_fit(
    geesolver_fit = geesolver_fit,
    xnames = xnames,
    qr_model_matrix = qr_model_matrix,
    family = family,
    weights = weights,
    y = y,
    model_matrix = model_matrix,
    model_frame = model_frame,
    id = id,
    repeated = repeated,
    call = call,
    data = data,
    model_terms = model_terms,
    control = control,
    method = method_original,
    association_structure = orstr
  )
  fit$converged <- (geesolver_fit$criterion[fit$iter] <= tolerance)
  if (method_original %in% c("bcgee-naive", "bcgee-robust", "bcgee-empirical")) {
    fit$converged <- TRUE
  }

  fit$alpha <- if (identical(orstr, "independence")) 0 else as.numeric(geesolver_fit$alpha)

  if (orstr %in% c("unstructured", "fixed")) {
    pairs_matrix <- combn(max(repeated), 2)
    alpha_names <- paste0("alpha_", pairs_matrix[1, ], ".", pairs_matrix[2, ])
    names(fit$alpha) <- alpha_names
  }
  if (!fit$converged) {
    warning("geewa_binary: algorithm did not converge", call. = FALSE)
  }
  eps <- 10 * .Machine$double.eps
  if (any(fit$fitted.values > 1 - eps) || any(fit$fitted.values < eps)) {
    warning("geewa_binary: fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
  }
  fit <- new_geer(fit)
  fit
}
