#' Solving Adjusted Generalized Estimating Equations
#'
#'
#' @description
#' Produces a Generalized Estimating Equations fit of the data and allows for
#' bias-reduction, bias-correction or penalized adjustment vectors.
#'
#'
#' @details
#' The \code{data} must be provided in case level or equivalently in `long'
#' format. See details about the `long' format in the function \link{reshape}.
#'
#' A term of the form \code{offset(expression)} is allowed in the right hand
#' side of \code{formula}.
#'
#' The default set for the \code{id} labels is \eqn{\{1,\ldots,N\}}, where
#' \eqn{N} is the sample size. If otherwise, the function recodes the given
#' labels of \code{id} onto this set.
#'
#' The argument \code{repeated} can be safely ignored only if the \code{data} is
#' written in such a way that the \eqn{t}-th observation in each cluster is
#' recorded at the \eqn{t}-th measurement occasion or when
#' \code{correlation_structure = "independence"}. Otherwise, the user
#' must provide \code{repeated}. The suggested set for the labels of
#' \code{repeated} is \eqn{\{1,\ldots,T\}}, where \eqn{T} is the number of
#' observed labels. If otherwise, the function recodes the given labels of
#' \code{repeated} onto this set.
#'
#' The variables \code{id} and \code{repeated} do not need to be pre-sorted.
#' Instead the function reshapes \code{data} in an ascending order of \code{id}
#' and \code{repeated}.
#'
#' The argument \code{method} specifies the adjustment added to the generalized
#' estimating equations. If \code{method = "gee"}, then the original GEE are
#' solved, i.e. no adjustment vector is added. If \code{method = "brgee_naive"},
#' \code{method = "brgee_robust"} or \code{method = "brgee_empirical"},
#' then the GEE will be adjusted to produce naive, robust or empirical
#' bias-reducing estimators, respectively. If \code{method = "bcgee_naive"},
#' \code{method = "bcgee_robust"} or \code{method = "bcgee_empirical"}, then the
#' GEE will be adjusted to produce naive, robust or emprical bias-corrected
#' estimators. If \code{method = "pgee_jeffreys"}, then penalized GEE will be
#' solved.
#'
#'
#' @param formula a formula expression of the form \code{response ~ predictors}
#'        as for other regression models. See the documentation of \link{glm} and
#'        \link{formula}.
#' @param data an optional data frame containing the variables provided in
#'        \code{formula}, \code{id} and \code{repeated}.
#' @param id a vector that identifies the clusters.
#' @param repeated an optional vector that identifies the order of observations
#'        within each cluster.
#' @param weights an optional vector of ‘prior weights’ to be used in the
#'        fitting process. Should be NULL or a numeric vector.
#' @param family a character, a family function of the results of a call to a
#'        family function describing the marginal distribution and the link function to
#'        be used in the marginal model.
#' @param control list that specifies the control variables for the GEE solver.
#' @param beta_start numerical vector indicating the initial values of the
#'        regression parameter vector.
#' @param correlation_structure a character indicating the working correlation
#'        structure. Options include \code{"independence"}, \code{"exchangeable"},
#'        \code{"ar1"}, \code{"m-dependent"} or \code{"unstructured"}.
#' @param use_p if set to \code{FALSE} then do not subtract the number of
#'        parameters when you are estimating phi or alpha.
#' @param Mv a positive integer which must be specified whenever
#'        \code{correlation_structure = "m-dependent"}.
#' @param alpha_vector numerical vector indicating the correlation structure for
#'        when \code{correlation_structure == "fixed"}. Otherwise, it is ignored.
#' @param phi_fixed logical indicating whether the phi value is fixed or not.
#' @param phi_value positive number indicating the value of phi when
#'        \code{phi_fixed == TRUE}.
#' @param method a character indicating the form of estimating equations.
#'        Options include \code{"gee"}, \code{"brgee_naive"},
#'        \code{"brgee_robust"}, \code{"brgee_empirical"}, \code{"bcgee_naive"},
#'        \code{"bcgee_robust"}, \code{"bcgee_empirical"} or \code{"pgee_jeffreys"}.
#' @param control_glm list of parameters for controlling the fitting process for the
#'        initial parameter vector values when these are not provided.
#' @param ... further arguments passed to/or from other methods.
#'
#'
#' @return \code{geewa} return an object from the class "geer". The function
#' \code{summary} (i.e., \code{summary.geer}) can be used to obtain or print a
#' summary of the results and the function \code{wald_test} to compare nested
#' models.
#'
#' The generic accessor functions \code{coefficients}, \code{fitted.values} and
#' \code{residuals} can be used to extract various useful features of the value
#' returned by \code{geer}.
#'
#' \item{call}{the matched call.}
#' \item{coefficients}{a named vector of coefficients.}
#' \item{naive_covariance}{the naive covariance matrix.}
#' \item{robust_covariance}{the robust covariance matrix.}
#' \item{bias_corrected_covariance}{the bias-corrected covariance matrix.}
#' \item{association_structure}{the name of the working assumption about the
#' association structure.}
#' \item{alpha}{a vector of the association parameters.}
#' \item{phi}{the scale parameter.}
#' \item{niter}{the number of iterations used.}
#' \item{criterion}{the value of the convergence criterion.}
#' \item{converged}{logical indicating whether the algorithm converged.}
#' \item{terms}{the \link{terms} object used.}
#' \item{y}{the response vector.}
#' \item{model_matrix}{the model matrix.}
#' \item{obs_no}{the number of observations across clusters.}
#' \item{family}{the \link{family} object used.}
#' \item{fitted.values}{a numeric vector of fitted values.}
#' \item{residuals}{a numeric vector of residuals.}
#' \item{linear.predictors}{a numeric vector of linear predictors.}
#' \item{clusters_no}{the number of clusters.}
#' \item{min_cluster_size}{the minimum cluster size.}
#' \item{max_cluster_size}{the maximum cluster size.}
#' \item{method}{the method that created the adjustment vector to the GEE.}
#'
#' @author Anestis Touloumis
#'
#' @examples
#' data("leprosy")
#' fitted_model_gee <- geewa(
#'    formula = bacilli ~ time+time:I(treatment == "C"),
#'    data = leprosy,
#'    id = id,
#'    family = poisson(link = "log"),
#'    correlation_structure = "exchangeable",
#'    method = "gee")
#' summary(fitted_model_gee, type = "bias-corrected")
#' fitted_model_brgee_robust <-
#'  update(fitted_model_gee, method = "brgee_robust")
#' summary(fitted_model_brgee_robust, type = "bias-corrected")
#' fitted_model_brgee_naive <-
#'  update(fitted_model_gee, method = "brgee_naive")
#' summary(fitted_model_brgee_naive, type = "bias-corrected")
#' fitted_model_brgee_empirical <-
#'   update(fitted_model_gee, method = "brgee_robust")
#' summary(fitted_model_brgee_empirical, type = "bias-corrected")#'
#'
#'
#' @export
geewa <-
  function(formula = formula(data),
           data = parent.frame(),
           id = id,
           repeated = NULL,
           family = gaussian(link = "identity"),
           weights,
           control = list(...),
           beta_start = NULL,
           correlation_structure = "independence",
           use_p = TRUE,
           Mv = 1,
           alpha_vector = NULL,
           phi_fixed = FALSE,
           phi_value = 1,
           method = "gee",
           control_glm = list(...),
           ...) {
    ## call
    if (missingArg(data))
      data <- environment(eval(formula))
    call <- match.call(expand.dots = TRUE)
    if (missing(data))
      data <- environment(formula)
    mcall <- match.call(expand.dots = FALSE)
    mnames <- c("formula", "data", "id", "repeated", "weights")
    mf <- match(mnames, names(mcall), 0L)
    m <- mcall[c(1L, mf)]
    mcall$drop.unused.levels <- TRUE
    m[[1L]] <- as.name("model.frame")
    model_frame <- eval(m, envir = parent.frame())

    if (is.null(m$id))
      m$id <- as.name("id")

    ## get explanatory variables
    model_terms <- attr(model_frame, "terms")

    ## extract response
    y <- model.response(model_frame)
    if (is.null(y))
      stop("response variable not found")
    y <- as.numeric(y)

    ## weights
    weights <- as.vector(model.weights(model_frame))
    if (is.null(weights)) {
      weights <- rep(1, length(y))
    } else {
      if (!is.numeric(weights))
        stop("'weights' must be a numeric vector")
      if (any(weights <= 0))
        stop("negative weights not allowed")
    }

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
      id_list <- split(id, id)
      order_id_index <- order(unlist(id_list))
      repeated <- unlist(lapply(id_list, function(x) seq.int(length(x))))
      repeated <- repeated[order_id_index]
    }
    repeated <- as.numeric(factor(repeated))

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
    offset <- model.extract(model_frame, "offset")
    if (length(offset) <= 1) offset <- rep(0, length(y))
    offset <- as.double(offset)

    ## check offset and y lengths
    if (length(offset) != length(y))
      stop("response variable and 'offset' are not of same length")

    ## extract model matrix
    model_matrix <- model.matrix(model_terms, model_frame)

    ## extract explanatory variable names
    xnames <- colnames(model_matrix)

    ## order response, model matrix, id and repeated
    ordered_index <- order(id, repeated)
    y <- y[ordered_index]
    model_matrix <- model_matrix[ordered_index, ]
    id <- id[ordered_index]
    repeated <- repeated[ordered_index]
    offset <- offset[ordered_index]

    ## convert model matrix to a matrix when p = 0
    if (length(xnames) == 1)
      model_matrix <- matrix(model_matrix, ncol = 1)

    ## check rank of model matrix
    qr_model_matrix <- qr(model_matrix)
    if (qr_model_matrix$rank < ncol(model_matrix))
      stop("rank-deficient model matrix")


    ## check correlation structure
    correlation_structures <- c("independence", "exchangeable", "ar1",
                                "m-dependent", "unstructured", "fixed")
    icheck <- as.integer(match(correlation_structure, correlation_structures,
                               -1))
    if (icheck < 1)
      stop("`correlation_structure` must be one of `independence`,
           `exchangeable`, `ar1`, `m-dependent`, `unstructured` or `fixed`")
    if (correlation_structure == "m-dependent" & ((Mv <= 0) | Mv %% 1 != 0))
      stop("Mv must be a positive integer number")
    if (correlation_structure !=  "m-dependent") Mv <- 1
    if (correlation_structure != "fixed") {
      alpha_vector <- 0
      alpha_fixed <- 0
    } else {
      if (is.null(alpha_vector))
        stop("`alpha_vector` must be provided when `correlation_structure == fixed`")
      alpha_vector <- as.numeric(alpha_vector)
      repeated_max <- max(repeated)
      if (length(alpha_vector) != choose(repeated_max, 2))
        stop("`alpha_vector` must be of a vector of size ",
             choose(repeated_max, 2))
      corr_matrix <- get_correlation_matrix(correlation_structure,
                                            alpha_vector,
                                            repeated_max)
      alpha_fixed <- 1
      }

    ## family
    familys <- c("gaussian", "poisson", "binomial", "Gamma", "inverse.gaussian",
                 "quasi", "quasibinomial", "quasipoisson")
    if (is.character(family))
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    icheck <- as.integer(match(family$family, familys, -1))
    if (icheck < 1)
      stop("`family` must be one of `gaussian`, `poisson`, `binomial`, `Gamma`,
           `inverse.gaussian`, `quasi`, `quasibinomial` or `quasipoisson`")

    if (family$family %in% c("quasi", "quasibinomial", "quasipoisson")) {
      if (family$family == "quasi") {
        family$family <- switch(family$varfun,
                                constant = "gaussian",
                                `mu(1-mu)` = "binomial",
                                mu = "poisson",
                                `mu^2` = "Gamma",
                                `mu^3` = "inverse.gaussian")
      } else {
        family$family <- switch(family$family,
                                quasibinomial = "binomial",
                                quasipoisson = "poisson")
      }
      family <- do.call(family$family, list(link = family$link))
    }

    ## link function
    links <- c("logit", "probit", "cauchit", "cloglog", "identity", "log",
               "sqrt", "1/mu^2", "inverse")
    link <- as.character(family$link)
    icheck <- as.integer(match(link, links, -1))
    if (icheck < 1)
      stop("`link` must be one of `logit`, `probit`, `cauchit`, `cloglog`,
           `identity, `log`, `sqrt``, `1/mu^2` or `inverse`")

    ## run fit
    methods <- c("gee",
                 "brgee_naive", "brgee_robust", "brgee_empirical",
                 "bcgee_naive", "bcgee_robust", "bcgee_empirical",
                 "pgee_jeffreys")
    method <- as.character(method)
    icheck <- as.integer(match(method, methods, -1))
    if (icheck < 1)
      stop("`method` must be one of `gee`, `brgee_naive`, `brgee_robust`,
           `brgee_empirical`, `bcgee_naive`, `bcgee_robust`, `bcgee_empirical`
           or `pgee_jeffreys`")


    ## control variable
    control <- do.call("geer_control", control)

    ## maxiter
    maxiter <- control$maxiter

    ## tolerance
    tolerance <- control$tolerance

    ## initial beta
    if (is.null(beta_start)) {
      control_glm <- do.call("brglm_control", control_glm)
      if (link != "identity") {
        if (method %in% c("gee",
                          "bcgee_naive", "bcgee_robust", "bcgee_empirical",
                          "brgee_robust", "brgee_empirical")) {
          type <- "ML"
        } else if (method  == "pgee_jeffreys") {
          type <- "MPL_Jeffreys"
        } else {
          type <- "AS_mean"
        }
        beta_zero <- brglmFit(x = model_matrix,
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
                                             a = control$jeffreys_power))$coefficients
      } else {
        beta_zero <- glm.fit(x = model_matrix,
                             y = y,
                             family = family,
                             weights = weights,
                             offset = offset,
                             control = list(epsilon = control$tolerance,
                                            maxit = control_glm$maxit,
                                            trace = FALSE))$coef
        }
      } else {
        beta_start <- as.numeric(beta_start)
        p <- ncol(model_matrix)
        if (length(beta_start) != p)
          stop("`beta_start` must be of a vector of length ", p)
        beta_zero <- beta_start
        }

    ## phi
    phi_fixed <- ifelse(phi_fixed, 1, 0)
    if (phi_fixed == 1) {
      phi_value <- as.numeric(phi_value)
      if (length(phi_value) < 0 & phi_value <= 0)
        stop("'phi_value' must be positive when `phi_fixed` is TRUE")
    } else {
      phi_value <- 1
    }



    ## use or not p when updating phi and alpha
    use_p <- as.logical(use_p)
    subtract_p <- ifelse(use_p, ncol(model_matrix), 0)

    method_original <- method
    if (method_original %in% c("bcgee_naive", "bcgee_robust", "bcgee_empirical")) {
      method <- "gee"
    }

    geesolver_fit <- fit_geesolver(y,
                                  model_matrix,
                                  id,
                                  repeated,
                                  link,
                                  family$family,
                                  beta_zero,
                                  correlation_structure,
                                  Mv,
                                  subtract_p,
                                  maxiter,
                                  tolerance,
                                  offset,
                                  phi_value,
                                  phi_fixed,
                                  alpha_vector,
                                  alpha_fixed,
                                  method,
                                  weights,
                                  control$step_maxiter,
                                  control$step_multiplier,
                                  control$jeffreys_power)

    if (method_original %in% c("bcgee_naive", "bcgee_robust", "bcgee_empirical")) {
      method <- sub("bcgee", "brgee", method_original)
      geesolver_fit <- fit_geesolver(y,
                                     model_matrix,
                                     id,
                                     repeated,
                                     link,
                                     family$family,
                                     as.numeric(c(geesolver_fit$beta_hat)),
                                     correlation_structure,
                                     Mv,
                                     subtract_p,
                                     1,
                                     tolerance,
                                     offset,
                                     phi_value,
                                     phi_fixed,
                                     alpha_vector,
                                     alpha_fixed,
                                     method,
                                     weights,
                                     1,
                                     1,
                                     control$jeffreys_power)
      method <- method_original
    }


    ## output
    fit <- list()

    fit$call <- call

    fit$coefficients <-  as.numeric(c(geesolver_fit$beta_hat))
    names(fit$coefficients) <- xnames

    fit$naive_covariance <- geesolver_fit$naive_covariance
    dimnames(fit$naive_covariance) <- list(xnames, xnames)
    fit$robust_covariance <- geesolver_fit$robust_covariance
    dimnames(fit$robust_covariance) <- list(xnames, xnames)
    fit$bias_corrected_covariance <- geesolver_fit$bc_covariance
    dimnames(fit$bias_corrected_covariance) <- list(xnames, xnames)

    fit$association_structure <- correlation_structure
    if (correlation_structure == "independence") {
      fit$alpha <- 0
      } else {
        fit$alpha <- c(geesolver_fit$alpha)
        }
    if (correlation_structure == "unstructured" | correlation_structure == "fixed") {
      pairs_matrix <- combn(max(repeated), 2)
      alpha_names <- paste0("alpha_",
                            paste(pairs_matrix[1, ], pairs_matrix[2, ], sep = ".")
                            )

      } else {
      alpha_names <- NULL
      }
    names(fit$alpha) <- alpha_names


    fit$phi <- geesolver_fit$phi



    fit$niter <- ncol(geesolver_fit$beta_mat) - 1
    fit$criterion <- geesolver_fit$criterion[fit$niter]
    fit$converged <- fit$criterion <= tolerance
    if (method_original %in% c("bcgee_naive", "bcgee_robust", "bcgee_empirical")) {
      fit$converged = TRUE
    }

    fit$terms <- model_terms
    fit$contrasts <- attr(model_matrix, "contrasts")
    fit$levels <- .getXlevels(attr(model_frame,"terms"), model_frame)
    fit$y <- y
    fit$model_matrix <- model_matrix
    fit$obs_no <- nrow(model_matrix)
    fit$family <- family

    fit$fitted.values <- c(geesolver_fit$fitted)
    fit$residuals <- c(geesolver_fit$residuals)
    fit$linear.predictors <- c(geesolver_fit$eta)

    fit$id <- as.numeric(id)
    fit$repeated <- as.numeric(repeated)
    fit$clusters_no <- max(id)
    clusters_sizes <- unlist(lapply(split(repeated, id), length))
    fit$min_cluster_size <- min(unique(clusters_sizes))
    fit$max_cluster_size <- max(unique(clusters_sizes))

    fit$method <- method
    fit$weights <- weights


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


    fit
  }

