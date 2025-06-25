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
#' respiratory2 <- respiratory[respiratory$center==2, ]
#' fitted_model <-
#'   geewa_binary(formula = y ~ baseline + treatment*gender + visit*age,
#'                id = id,
#'                repeated = visit,
#'                link = "probit",
#'                data = respiratory2,
#'                orstr = "independence",
#'                method = "pgee-jeffreys")
#' summary(fitted_model, type = "bias-corrected")
#'
#' data("obesity")
#' fitted_model_gee <- geewa_binary(
#'   formula = I(obesity == "Obese") ~ age + zygosity + ancestry + bacteroides,
#'   id = fid,
#'   data = obesity,
#'   link = "logit",
#'   orstr = "independence",
#'   method = "gee")
#' summary(fitted_model_gee, type = "bias-corrected")
#' fitted_model_brgee_robust <-
#'   update(fitted_model_gee, method = "brgee-robust")
#' summary(fitted_model_brgee_robust, type = "bias-corrected")
#' fitted_model_brgee_naive <-
#'   update(fitted_model_gee, method = "brgee-naive")
#' summary(fitted_model_brgee_naive, type = "bias-corrected")
#' fitted_model_brgee_empirical <-
#'   update(fitted_model_gee, method = "brgee-robust")
#' summary(fitted_model_brgee_empirical, type = "bias-corrected")
#'
#' @export
geewa_binary <- function(formula = formula(data),
                         link = "logit",
                         data = parent.frame(),
                         id = id,
                         repeated = NULL,
                         control = geer_control(),
                         orstr = "exchangeable",
                         method = "gee",
                         weights,
                         beta_start = NULL,
                         control_glm = list(...),
                         alpha_vector = NULL,
                         ...) {
  ## call
  call <- match.call(expand.dots = TRUE)
  if (missing(data))
    data <- environment(formula)
  mcall <- match.call(expand.dots = FALSE)
  mcall$drop.unused.levels <- TRUE
  mnames <- c("formula", "data", "id", "repeated", "weights", "offset")
  mf <- match(mnames, names(mcall), 0L)
  m <- mcall[c(1L, mf)]
  m[[1]] <- as.name("model.frame")
  model_frame <- eval(m, envir = parent.frame())
  ## link function
  links <- c("logit", "probit", "cauchit", "cloglog", "identity", "log",
             "sqrt", "1/mu^2", "inverse")
  link <- as.character(link)
  icheck <- as.integer(match(link, links, -1))
  if (icheck < 1)
    stop("'link' should be one of 'logit', 'probit', 'cauchit', 'cloglog',
           'identity', 'log', 'sqrt', '1/mu^2', 'inverse'")
  family <- binomial(link = link)
  ## extract response
  y <- model.response(model_frame, "any")
  if (is.null(y))
    stop("response variable not found")
  if (any(class(y) %in% c("factor", "character"))) {
    y <- as.numeric(y != levels(factor(y))[1])
  }
  ## weights
  weights <- as.vector(model.weights(model_frame))
  if (is.null(weights)) {
    weights <- rep(1, length(y))
  } else {
    if (!is.numeric(weights))
      stop("'weights' should be a numeric vector")
    if (any(weights <= 0))
      stop("non-positive weights not allowed")
  }
  if (is.matrix(y) && ncol(y) == 2) {
    weights <- apply(y, 1, sum)
    y <- y[, 1]/weights
  }
  y <- as.numeric(y)
  ## extract id and map values to 1,.., N
  id <- model.extract(model_frame, "id")
  if (is.null(id))
    stop("'id' not found")
  id <- as.numeric(factor(id))
  ## check id and y lengths
  if (length(id) != length(y))
    stop("response variable and 'id' are not of same length")
  ## extract or create repeated and map values to 1,..., T
  repeated <- model.extract(model_frame, "repeated")
  if (is.null(repeated)) {
    id_list <- split(id, id)[unique(id)]
    repeated <- unlist(lapply(id_list, function(x) seq.int(length(x))))
  } else {
    repeated <- as.numeric(factor(repeated))
  }
  ## check repeated and y lengths
  if (length(repeated) != length(y))
    stop("response variable and 'repeated' are not of same length")
  ## check if id and repeated are identical
  if (all(id == repeated))
    stop("'repeated' and 'id' must be different")
  ## check duplicated values in repeated for each subject
  if (any(unlist(lapply(split(repeated, id), duplicated))))
    stop("'repeated' does not have unique values per 'id'")
  ## offset term
  offset <- model.offset(model_frame)
  if (length(offset) <= 1) offset <- rep(0, length(y))
  offset <- as.double(offset)
  ## check offset and y lengths
  if (length(offset) != length(y))
    stop("response variable and 'offset' are not of same length")
  ## get model terms
  model_terms <- attr(model_frame, "terms")
  ## extract model matrix
  model_matrix <- model.matrix(model_terms, model_frame)
  ## extract explanatory variable names
  xnames <- colnames(model_matrix)
  ## convert model matrix to a matrix when p = 1
  if (length(xnames) == 1)
    model_matrix <- matrix(model_matrix, ncol = 1)
  ## check rank of model matrix
  qr_model_matrix <- qr(model_matrix)
  if (qr_model_matrix$rank < ncol(model_matrix))
    stop("rank-deficient model matrix")
  ## control variables
  control <- do.call("geer_control", control)
  maxiter <- control$maxiter
  tolerance <- control$tolerance
  ## estimation method
  methods <- c("gee",
               "brgee-naive", "brgee-robust", "brgee-robust2", "brgee-empirical",
               "bcgee-naive", "bcgee-robust", "bcgee-robust2", "bcgee-empirical",
               "pgee-jeffreys")
  method <- as.character(method)
  icheck <- as.integer(match(method, methods, -1))
  if (icheck < 1)
    stop("'method' should be one of 'gee', 'brgee-naive', 'brgee-robust', 'brgee-empirical', 'bcgee-naive', 'bcgee-robust', 'bcgee-empirical', 'pgee-jeffreys'")
  ## initial beta
  if (is.null(beta_start)) {
    control_glm <- do.call("brglm_control", control_glm)
    if (link != "identity") {
      if (method %in% c("gee", "bcgee-robust2",
                        "bcgee-naive", "bcgee-robust", "bcgee-empirical",
                        "brgee-robust", "brgee-empirical")) {
        type <- "ML"
      } else if (method  == "pgee-jeffreys") {
        type <- "MPL_Jeffreys"
      } else {
        type <- "AS_mean"
      }
      glmfit <- try(brglmFit(x = model_matrix,
                             y = y,
                             family = family,
                             weights = weights,
                             offset = offset,
                             control = list(epsilon = control$tolerance,
                                            maxit = control_glm$maxit,
                                            type = type,
                                            trace = FALSE,
                                            slowit = control_glm$slowit,
                                            max_step_factor = control_glm$max_step_factor,
                                            a = control$jeffreys_power)),
                    silent = TRUE)
    } else {
      glmfit <- try(glm.fit(x = model_matrix,
                            y = y,
                            family = family,
                            weights = weights,
                            offset = offset,
                            control = list(epsilon = control$tolerance,
                                           maxit = control_glm$maxit,
                                           trace = FALSE)),
                    silent = TRUE)
    }
    if (!inherits(glmfit, "try-error")) {
      beta_zero <- glmfit$coefficients
    } else {
      stop("cannot find valid starting values: please specify some!!",
           call. = FALSE)
    }
  } else {
    beta_start <- as.numeric(beta_start)
    p <- ncol(model_matrix)
    if (length(beta_start) != p)
      stop("'beta_start' must be of a vector of length ", p)
    beta_zero <- beta_start
  }
  ## estimating marginalized odds ratios
  orstrs <- c("independence", "exchangeable", "unstructured", "fixed")
  icheck <- as.integer(match(orstr, orstrs, -1))
  if (icheck < 1)
    stop("'orstr' should be one of 'independence', 'exchangeable', 'unstructured', 'fixed'")
  if (orstr == "independence") {
    adding_constant <- NULL
    alpha_vector <- rep.int(1, choose(max(repeated), 2))
  } else if (orstr == "fixed") {
    adding_constant <- NULL
    pairs_no <- choose(max(repeated), 2)
    if (!is.numeric(alpha_vector))
      stop("'alpha_vector' must be as a numeric vector")
    if (length(alpha_vector) != pairs_no)
      stop("'alpha_vector' must be of length ", pairs_no)
  } else {
    adding_constant <- control$or_adding
    alpha_vector <- get_marginalized_odds_ratios(round(y),
                                                 id,
                                                 repeated,
                                                 weights,
                                                 adding_constant,
                                                 orstr)
  }
  ## change the estimation method to gee for bias corrected estimators
  method_original <- method
  if (method_original %in% c("bcgee-naive", "bcgee-robust", "bcgee-robust2", "bcgee-empirical")) {
    method <- "gee"
  }
  ## gee with or without adjustments
  geesolver_fit <- fit_bingee_or(y, model_matrix, id, repeated, weights, link,
                                 beta_zero, offset, maxiter, tolerance,
                                 control$step_maxit, control$step_multi,
                                 control$jeffreys_power, method, alpha_vector)
  ## only for bias-corrected estimators
  if (method_original %in% c("bcgee-naive", "bcgee-robust", "bcgee-robust2", "bcgee-empirical")) {
    if (geesolver_fit$criterion[ncol(geesolver_fit$beta_mat) - 1] <= tolerance) {
      method <- sub("bcgee", "brgee", method_original)
      geesolver_fit <- fit_bingee_or(y, model_matrix, id, repeated, weights,
                                     link, as.numeric(c(geesolver_fit$beta_hat)),
                                     offset, 1, tolerance, 1, 1,
                                     control$jeffreys_power, method, alpha_vector)
      method <- method_original
    } else {
      stop("bias-corrected GEE estimator is NA due to non-convergence of the GEE model")
    }
  }
  ## output
  fit <- list()
  fit$coefficients <-  as.numeric(c(geesolver_fit$beta_hat))
  names(fit$coefficients) <- xnames
  fit$residuals <- c(geesolver_fit$residuals)
  fit$fitted.values <- c(geesolver_fit$fitted)
  fit$rank <- qr(model_matrix)$rank
  fit$family <- family
  fit$linear.predictors <- c(geesolver_fit$eta)
  fit$iter <- ncol(geesolver_fit$beta_mat) - 1
  fit$prior.weights <- weights
  fit$df.residual <- nrow(model_matrix) - ncol(model_matrix)
  fit$y <- y
  fit$x <- model_matrix
  fit$id <- as.numeric(id)
  fit$repeated <- as.numeric(repeated)
  fit$converged <- geesolver_fit$criterion[fit$iter] <= tolerance
  if (method_original %in% c("bcgee-naive", "bcgee-robust", "bcgee-robust2", "bcgee-empirical")) {
    fit$converged <- TRUE
  }
  fit$call <- call
  fit$formula <- fit$call$formula
  fit$terms <- model_terms
  fit$data <- data
  fit$offset <- geesolver_fit$offset
  fit$control <- control
  fit$method <- method
  fit$contrasts <- attr(model_matrix, "contrasts")
  fit$xlevels <- .getXlevels(attr(model_frame,"terms"), model_frame)
  fit$naive_covariance <- geesolver_fit$naive_covariance
  dimnames(fit$naive_covariance) <- list(xnames, xnames)
  fit$robust_covariance <- geesolver_fit$robust_covariance
  dimnames(fit$robust_covariance) <- list(xnames, xnames)
  fit$bias_corrected_covariance <- geesolver_fit$bc_covariance
  dimnames(fit$bias_corrected_covariance) <- list(xnames, xnames)
  fit$association_structure <- orstr
  if (orstr == "independence") {
    fit$alpha <- 1
  } else if (orstr == "exchangeable") {
    fit$alpha <- c(geesolver_fit$alpha)[1]
  } else {
    fit$alpha <- c(geesolver_fit$alpha)
  }
  if (orstr == "unstructured" | orstr == "fixed") {
    pairs_matrix <- combn(max(repeated), 2)
    alpha_names <- paste0("alpha_",
                          paste(pairs_matrix[1, ], pairs_matrix[2, ], sep = ".")
    )
  } else {
    alpha_names <- NULL
  }
  names(fit$alpha) <- alpha_names
  fit$phi <- 1
  fit$obs_no <- nrow(model_matrix)
  fit$clusters_no <- max(id)
  clusters_sizes <- unlist(lapply(split(repeated, id), length))
  fit$min_cluster_size <- min(unique(clusters_sizes))
  fit$max_cluster_size <- max(unique(clusters_sizes))
  class(fit) <- "geer"
  if (!fit$converged)
    warning("geewa: algorithm did not converged",
            call. = FALSE)
  eps <- 10 * .Machine$double.eps
  if (any(fit$fitted.values > 1 - eps) || any(fit$fitted.values < eps))
    warning("geewa_binary: fitted probabilities numerically 0 or 1 occurred",
            call. = FALSE)
  fit
}
