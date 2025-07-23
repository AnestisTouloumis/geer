#' @title
#' Solving (Adjusted) GEE
#'
#' @description
#' Produces a(n) (adjusted) Generalized Estimating Equations (GEE) fit of the
#' data. The estimation methods include the traditional GEE method,
#' bias-reduction methods, bias-correction methods and a penalized GEE method.
#'
#' @param formula \code{formula} expression of the form
#'        \code{response ~ predictors}: a symbolic description of the marginal
#'        model to be fitted.
#' @param family \code{\link{family}} object indicating the link and variance
#'        functions. Options include \code{gaussian}, \code{binomial},
#'        \code{poisson}, \code{Gamma}, \code{inverse.gaussian} and \code{quasi}.
#'         By default, \code{family = gaussian(link = "identity")}.
#' @param data optional data frame containing the variables provided in
#'        \code{formula}, \code{id} and \code{repeated}.
#' @param id a vector identifying the clusters.
#' @param repeated optional vector identifying the order of observations
#'        within each cluster.
#' @param control list specifying the control variables for the GEE solver
#'        and the power of the Jeffreys prior when
#'        \code{method = "pgee-jeffreys"}.
#' @param corstr character specifying the correlation structure. Options
#'        include \code{"independence"}, \code{"exchangeable"}, \code{"ar1"},
#'        \code{"m-dependent"} and \code{"unstructured"}. By default,
#'        \code{corstr = "independence"}.
#' @param Mv positive integer which must be specified whenever
#'        \code{corstr = "m-dependent"}. By default, \code{Mv = 1}.
#' @param method character specifying the (adjusted) GEE. Options include the
#'        traditional gee method (\code{"gee"}), bias-reduction methods
#'        (\code{"brgee-naive"}, \code{"brgee-robust"},
#'        \code{"brgee-empirical"}), bias-correction methods
#'        (\code{"bcgee-naive"}, \code{"bcgee-robust"},
#'        \code{"bcgee-empirical"}) and a penalized GEE method
#'        (\code{"pgee-jeffreys"}). By default, \code{method = "gee"}.
#' @param weights optional numeric vector identifying the weights for each
#'        observation.
#' @param beta_start numerical vector indicating the initial values of the
#'        regression parameters.
#' @param control_glm list of parameters for controlling the fitting process for
#'        the initial values of the parameter vector when these are not provided.
#' @param use_p logical indicating whether to use the \code{N - p} correction
#'        for estimating the scale and the correlation parameters. By default,
#'        \code{use_p = TRUE}.
#' @param alpha_vector numerical vector indicating the correlation structure for
#'        when \code{corstr == "fixed"}. Otherwise, it is ignored.
#' @param phi_fixed logical indicating whether the scale parameter is fixed at
#'        the value of \code{phi_value}. By default, \code{phi_fixed = FALSE}.
#' @param phi_value positive number indicating the value to which the scale
#'        parameter must be fixed; only used when \code{phi_fixed == TRUE}.
#'        By default, \code{phi_value = 1}.
#' @param ... further arguments passed to/or from other methods.
#'
#' @details
#' \code{method} specifies the adjustment vector (if any) added to the ordinary
#' generalized estimating equations. If \code{method = "gee"}, then the original
#' GEE are solved, i.e. no adjustment vector is added. If
#' \code{method = "brgee-naive"}, \code{method = "brgee-robust"} or
#' \code{method = "brgee-empirical"}, then the GEE will be adjusted to
#' produce naive, robust or empirical bias-reducing estimators, respectively.
#' If \code{method = "bcgee-naive"}, \code{method = "bcgee-robust"} or
#' \code{method = "bcgee-empirical"}, then the GEE will be adjusted to
#' produce naive, robust or empirical bias-corrected estimators. If
#' \code{method = "pgee-jeffreys"}, then the GEE will be penalized using a
#' Jeffreys prior type penalty.
#'
#' For the construction of the \code{formula} argument, see the documentation of
#' \code{\link{glm}} and \code{\link{formula}}.
#'
#' The \code{data} must be provided in case level or equivalently in `long'
#' format. See details about the `long' format in the function
#' \code{\link{reshape}}.
#'
#' The default set for the \code{id} labels is \eqn{\{1,\ldots,N\}}, where
#' \eqn{N} is the sample size, i.e. the total number of clusters. Otherwise,
#' the function recodes the given labels of \code{id} onto this set.
#'
#' The argument \code{repeated} can be safely ignored only if the \code{data}
#' are assumed to be sorted so that observations on a cluster are contiguous
#' rows for all entities in the formula. If this is not the case, then the user
#' is advised to provide \code{repeated}. The default set for the labels of
#' \code{repeated} is \eqn{\{1,\ldots,T\}}, where \eqn{T} is the number of
#' observed labels. Otherwise, the function recodes the given labels of
#' \code{id} onto this set.
#'
#' The variables \code{id} and \code{repeated} do not need to be pre-sorted.
#' Instead the function reshapes \code{data} in an ascending order of \code{id}
#' and \code{repeated}.
#'
#' A term of the form \code{offset(expression)} is allowed in the right hand
#' side of \code{formula}.
#'
#' The length of \code{id} and, of \code{repeated} or \code{weight} when these
#' are provided, should be the same as the number of observations.
#'
#' @return
#' Returns an object from the class \code{geer}, a list with components:
#' \item{coefficients}{a named vector of coefficients.}
#' \item{residuals}{the working residuals.}
#' \item{fitted.values}{the fitted mean values, obtained by transforming the
#'      linear predictors by the inverse of the link function.}
#' \item{rank}{the numeric rank of the fitted model.}
#' \item{family}{the \code{\link{family}} object used.}
#' \item{linear.predictors}{the linear fit on link scale.}
#' \item{iter}{the number of iterations used.}
#' \item{prior.weights}{the weights initially supplied, a vector of 1s if none
#'       were.}
#' \item{df.residual}{the residual degrees of freedom.}
#' \item{y}{the response vector.}
#' \item{x}{the model matrix.}
#' \item{id}{the \code{id} vector.}
#' \item{repeated}{the \code{repeated} vector.}
#' \item{converged}{logical indicating whether the algorithm converged.}
#' \item{call}{the matched call.}
#' \item{formula}{the formula supplied.}
#' \item{terms}{the \code{\link{terms}} object used.}
#' \item{data}{the data argument.}
#' \item{offset}{the offset vector used.}
#' \item{control}{the value of the control argument used.}
#' \item{method}{the method that created the adjustment vector to the GEE, if
#'       any.}
#' \item{contrasts}{the contrasts used.}
#' \item{xlevels}{a record of the levels of the factors used in fitting.}
#' \item{naive_covariance}{the naive covariance matrix.}
#' \item{robust_covariance}{the robust covariance matrix.}
#' \item{bias_corrected_covariance}{the bias-corrected covariance matrix.}
#' \item{association_structure}{the name of the working assumption about the
#'       association structure.}
#' \item{alpha}{a vector of the association parameters.}
#' \item{phi}{the scale parameter.}
#' \item{obs_no}{the number of observations across clusters.}
#' \item{clusters_no}{the number of clusters.}
#' \item{min_cluster_size}{the minimum cluster size.}
#' \item{max_cluster_size}{the maximum cluster size.}
#'
#' @author Anestis Touloumis
#'
#' @seealso \code{link{geewa_binary}}.
#'
#' @examples
#' data("leprosy")
#' fitted_model_gee <- geewa(formula = seizures ~ treatment + log(base/4) + log(age),
#'                           data = epilepsy, id = id, family = poisson(link = "log"),
#'                           corstr = "exchangeable", method = "gee")
#' summary(fitted_model_gee, type = "bias-corrected")
#' fitted_model_brgee_robust <- update(fitted_model_gee, method = "brgee-robust")
#' summary(fitted_model_brgee_robust, type = "bias-corrected")
#' fitted_model_brgee_naive <- update(fitted_model_gee, method = "brgee-naive")
#' summary(fitted_model_brgee_naive, type = "bias-corrected")
#' fitted_model_brgee_empirical <- update(fitted_model_gee, method = "brgee-robust")
#' summary(fitted_model_brgee_empirical, type = "bias-corrected")
#'
#' @export
geewa <- function(formula = formula(data),
                  family = gaussian(link = "identity"),
                  data = parent.frame(),
                  id = id,
                  repeated,
                  control = geer_control(...),
                  corstr = "independence",
                  Mv = 1,
                  method = "gee",
                  weights,
                  beta_start = NULL,
                  control_glm = list(...),
                  use_p = TRUE,
                  alpha_vector = NULL,
                  phi_fixed = FALSE,
                  phi_value = 1,
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
  mcall$drop.unused.levels <- TRUE
  m[[1L]] <- as.name("model.frame")
  model_frame <- eval(m, envir = parent.frame())
  ## family
  familys <- c("gaussian", "poisson", "binomial", "Gamma", "inverse.gaussian",
               "quasi")
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  icheck <- as.integer(match(family$family, familys, -1))
  if (icheck < 1)
    stop("'family' should be one of 'gaussian', 'poisson', 'binomial', 'Gamma',
           'inverse.gaussian', 'quasi'")
  links <- c("logit", "probit", "cauchit", "cloglog", "identity", "log",
             "sqrt", "1/mu^2", "inverse")
  link <- as.character(family$link)
  icheck <- as.integer(match(link, links, -1))
  if (icheck < 1)
    stop("'link' should be one of 'logit', 'probit', 'cauchit', 'cloglog',
           'identity', 'log', 'sqrt', '1/mu^2', 'inverse'")
  ## extract response
  y <- model.response(model_frame, "any")
  if (is.null(y))
    stop("response variable not found")
  if (family$family == "binomial" & any(class(y) %in% c("factor", "character"))) {
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
  if (family$family == "binomial") {
    if (is.matrix(y) && ncol(y) == 2) {
      weights <- apply(y, 1, sum)
      y <- y[, 1]/weights
    }
  }
  y <- as.numeric(y)
  ## extract id and map labels to 1,.., N
  id <- model.extract(model_frame, "id")
  if (is.null(id))
    stop("'id' not found")
  id <- as.numeric(factor(id))
  ## check id and y lengths
  if (length(id) != length(y))
    stop("response variable and 'id' are not of same length")
  ## extract or create repeated and map labels to 1,..., T
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
  ## get explanatory variables
  model_terms <- attr(model_frame, "terms")
  ## extract model matrix
  model_matrix <- model.matrix(model_terms, model_frame)
  ## extract explanatory variable names
  xnames <- colnames(model_matrix)
  ## convert model matrix to a matrix when p = 0
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
    stop("'method' should be one of 'gee',
               'brgee-naive', 'brgee-robust', 'brgee-empirical',
               'bcgee-naive', 'bcgee-robust', 'bcgee-empirical',
               'pgee-jeffreys'")
  ## initial beta
  if (is.null(beta_start)) {
    control_glm <- do.call("brglm_control", control_glm)
    if (link != "identity" & family$family != "quasi") {
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
      stop("'beta_start' should be of a vector of length ", p)
    beta_zero <- beta_start
  }
  ## check phi
  phi_fixed <- ifelse(phi_fixed, 1, 0)
  if (phi_fixed == 1) {
    phi_value <- as.numeric(phi_value)
    if (length(phi_value) < 0 & phi_value <= 0)
      stop("'phi_value' should be positive when `phi_fixed` is TRUE")
  } else {
    phi_value <- 1
  }
  ## check correlation structure
  corstrs <- c("independence", "exchangeable", "ar1",
               "m-dependent", "unstructured", "fixed")
  icheck <- as.integer(match(corstr, corstrs, -1))
  if (icheck < 1)
    stop("'corstr' should be one of 'independence',
           'exchangeable', 'ar1', 'm-dependent', 'unstructured', 'fixed'")
  if (corstr !=  "m-dependent") Mv <- 1
  if (corstr == "m-dependent" & ((Mv <= 0) | Mv %% 1 != 0))
    stop("Mv must be a positive integer number")
  if (corstr != "fixed") {
    alpha_vector <- 0
    alpha_fixed <- 0
  } else {
    if (is.null(alpha_vector))
      stop("'alpha_vector' should be provided when 'corstr == fixed'")
    alpha_vector <- as.numeric(alpha_vector)
    repeated_max <- max(repeated)
    if (length(alpha_vector) != choose(repeated_max, 2))
      stop("'alpha_vector' must be of a vector of size ", choose(repeated_max, 2))
    corr_matrix <- get_correlation_matrix(corstr,
                                          alpha_vector,
                                          repeated_max)
    alpha_fixed <- 1
  }
  ## use p in the denominator when updating phi and alpha
  use_p <- as.logical(use_p)
  subtract_p <- ifelse(use_p, ncol(model_matrix), 0)
  ## change the estimation method to gee for bias corrected estimators
  method_original <- method
  if (method_original %in% c("bcgee-naive", "bcgee-robust", "bcgee-robust2", "bcgee-empirical")) {
    method <- "gee"
  }

  if (family$family == "quasi") {
    family$family <- switch(family$varfun,
                            constant = "gaussian",
                            `mu(1-mu)` = "binomial",
                            mu = "poisson",
                            `mu^2` = "Gamma",
                            `mu^3` = "inverse.gaussian")
    family <- do.call(family$family, list(link = family$link))
    link <- family$link
  }


  ## gee with or without adjustments
  geesolver_fit <- fit_geesolver_cc(y, model_matrix, id, repeated, weights,
                                    link, family$family, beta_zero, offset,
                                    maxiter, tolerance, control$step_maxit,
                                    control$step_multi, control$jeffreys_power,
                                    method, subtract_p, alpha_vector, alpha_fixed,
                                    corstr, Mv, phi_value, phi_fixed)
  ## only for bias-corrected estimators
  if (method_original %in% c("bcgee-naive", "bcgee-robust", "bcgee-robust2", "bcgee-empirical")) {
    if (geesolver_fit$criterion[ncol(geesolver_fit$beta_mat) - 1] <= tolerance) {
      method <- sub("bcgee", "brgee", method_original)
      geesolver_fit <- fit_geesolver_cc(y, model_matrix, id, repeated, weights,
                                        link, family$family,
                                        as.numeric(c(geesolver_fit$beta_hat)),
                                        offset, 1, tolerance, 1, 1,
                                        control$jeffreys_power, method,
                                        subtract_p, geesolver_fit$alpha, 1,
                                        corstr, Mv, geesolver_fit$phi, 1)
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
  fit$association_structure <- corstr
  if (corstr == "independence") {
    fit$alpha <- 0
  } else {
    fit$alpha <- c(geesolver_fit$alpha)
  }
  if (corstr == "unstructured" | corstr == "fixed") {
    pairs_matrix <- combn(max(repeated), 2)
    alpha_names <- paste0("alpha_",
                          paste(pairs_matrix[1, ], pairs_matrix[2, ], sep = ".")
    ) } else {
      alpha_names <- NULL
    }
  names(fit$alpha) <- alpha_names
  fit$phi <- geesolver_fit$phi
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
  if (family$family == "binomial") {
    if (any(fit$fitted.values > 1 - eps) || any(fit$fitted.values < eps))
      warning("geewa: fitted probabilities numerically 0 or 1 occurred",
              call. = FALSE)
  }
  if (family$family == "poisson") {
    if (any(fit$fitted.values < eps))
      warning("geewa: fitted rates numerically 0 occurred",
              call. = FALSE)
  }
  if (any(eigen(get_correlation_matrix(corstr,fit$alpha, max(repeated)),
                symmetric = TRUE,
                only.values = TRUE)$values <= 0))
    warning("geewa: working correlation matrix is positive definite",
            call. = FALSE)
  fit
}
