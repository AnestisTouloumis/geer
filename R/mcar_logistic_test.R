make_mcar_covariate_formula <- function(object, formula) {
  if (is.null(formula)) {
    object_formula <- stats::formula(object)
    term_labels <- attr(object$terms, "term.labels")
    if (length(term_labels)) {
      formula <- stats::reformulate(term_labels)
    } else {
      formula <- stats::as.formula("~ 1")
    }
    environment(formula) <- environment(object_formula)
  } else {
    if (!inherits(formula, "formula")) {
      stop("'formula' must be a one-sided formula", call. = FALSE)
    }
    if (length(formula) != 2L) {
      stop(
        "'formula' must be one-sided, for example '~ treatment + age'",
        call. = FALSE
      )
    }
  }

  trms <- stats::terms(formula)
  if (!identical(attr(trms, "intercept"), 1L)) {
    stop(
      "'formula' must include an intercept; use '~ 1' when no covariates are required",
      call. = FALSE
    )
  }
  if (!is.null(attr(trms, "offset")) && length(attr(trms, "offset"))) {
    stop("offset terms are not supported in 'formula'", call. = FALSE)
  }

  formula
}


normalize_mcar_response <- function(response, object) {
  if (is.matrix(response) || is.data.frame(response)) {
    stop(
      "the regression-based MCAR diagnostic currently requires a univariate response at each occasion",
      call. = FALSE
    )
  }

  is_binomial_like <-
    identical(object$family$family, "binomial") ||
    identical(object$family$family, "quasibinomial") ||
    (identical(object$family$family, "quasi") &&
       !is.null(object$family$varfun) &&
       identical(object$family$varfun, "mu(1-mu)"))

  if (is_binomial_like && is.factor(response)) {
    if (nlevels(response) != 2L) {
      stop("a binomial factor response must have exactly two levels", call. = FALSE)
    }
    out <- as.numeric(response == levels(response)[2L])
    out[is.na(response)] <- NA_real_
    return(out)
  }
  if (is_binomial_like && is.character(response)) {
    observed <- response[!is.na(response)]
    lev <- levels(factor(observed))
    if (length(lev) != 2L) {
      stop("a binomial character response must have exactly two levels", call. = FALSE)
    }
    out <- as.numeric(response == lev[2L])
    out[is.na(response)] <- NA_real_
    return(out)
  }

  if (!is.numeric(response) && !is.logical(response)) {
    stop(
      "the response must be numeric or a two-level binary response for the regression-based MCAR diagnostic",
      call. = FALSE
    )
  }
  out <- as.numeric(response)
  if (any(!is.na(out) & !is.finite(out))) {
    stop("observed response values must be finite", call. = FALSE)
  }
  out
}


reconstruct_mcar_transition_frame <- function(object, formula, data = NULL) {
  if (is.null(data)) data <- object$data
  if (is.null(data)) {
    stop(
      "the original data are not available in the fitted object; supply them through 'data'",
      call. = FALSE
    )
  }

  object_formula <- stats::formula(object)
  formula_env <- environment(object_formula)
  if (is.null(formula_env)) formula_env <- parent.frame()
  eval_env <- new.env(parent = formula_env)
  assign(".geer_mcar_data", data, envir = eval_env)
  assign(".geer_original_formula", object_formula, envir = eval_env)

  response_call <- object$call
  response_call <- response_call[c(1L, match("data", names(response_call), 0L))]
  response_call[[1L]] <- quote(stats::model.frame)
  response_call$formula <- quote(.geer_original_formula)
  response_call$data <- quote(.geer_mcar_data)
  response_call$na.action <- quote(stats::na.pass)
  response_call$drop.unused.levels <- FALSE

  response_frame <- tryCatch(
    eval(response_call, envir = eval_env, enclos = formula_env),
    error = function(e) {
      stop(
        sprintf(
          "could not reconstruct the original response for the regression-based MCAR diagnostic: %s",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
  response <- normalize_mcar_response(
    stats::model.response(response_frame),
    object
  )
  if (!length(response)) {
    stop("the original response could not be recovered", call. = FALSE)
  }
  assign(".geer_mcar_response", response, envir = eval_env)

  rhs <- formula[[2L]]
  analysis_formula <- as.call(
    list(as.name("~"), as.name(".geer_mcar_response"), rhs)
  )
  class(analysis_formula) <- "formula"
  environment(analysis_formula) <- eval_env
  assign(".geer_mcar_formula", analysis_formula, envir = eval_env)

  mcall <- object$call
  keep <- c("data", "subset", "id", "repeated")
  mcall <- mcall[c(1L, match(keep, names(mcall), 0L))]
  mcall[[1L]] <- quote(stats::model.frame)
  mcall$formula <- quote(.geer_mcar_formula)
  mcall$data <- quote(.geer_mcar_data)
  mcall$na.action <- quote(stats::na.pass)
  mcall$drop.unused.levels <- TRUE

  model_frame <- tryCatch(
    eval(mcall, envir = eval_env, enclos = formula_env),
    error = function(e) {
      stop(
        sprintf(
          "could not construct the risk-set data for the regression-based MCAR diagnostic: %s",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  y <- stats::model.response(model_frame)
  if (length(y) != nrow(model_frame)) {
    stop("the response could not be reconstructed correctly", call. = FALSE)
  }

  id_raw <- stats::model.extract(model_frame, "id")
  if (is.null(id_raw)) {
    stop("'id' could not be recovered from the fitted model", call. = FALSE)
  }
  if (anyNA(id_raw)) {
    stop("'id' cannot contain missing values in the missingness model", call. = FALSE)
  }
  id <- as.numeric(factor(id_raw))

  repeated_raw <- stats::model.extract(model_frame, "repeated")
  if (is.null(repeated_raw)) {
    repeated <- stats::ave(id, id, FUN = seq_along)
  } else {
    if (anyNA(repeated_raw)) {
      stop("'repeated' cannot contain missing values in the missingness model", call. = FALSE)
    }
    repeated <- as.numeric(factor(repeated_raw))
  }

  if (length(y) != length(id) || length(y) != length(repeated)) {
    stop("response, 'id', and 'repeated' must have the same length", call. = FALSE)
  }
  if (any(duplicated(data.frame(id = id, repeated = repeated)))) {
    stop("'repeated' must identify unique measurements within each 'id'", call. = FALSE)
  }

  trms <- stats::terms(analysis_formula, data = model_frame)
  design <- stats::model.matrix(trms, model_frame)
  intercept <- which(colnames(design) == "(Intercept)")
  if (length(intercept) != 1L) {
    stop("the covariate formula must contain exactly one intercept", call. = FALSE)
  }
  covariate_index <- setdiff(seq_len(ncol(design)), intercept)
  covariates <- design[, covariate_index, drop = FALSE]
  if (anyNA(covariates) || any(!is.finite(covariates))) {
    stop(
      "covariates in the missingness model must be fully observed and finite",
      call. = FALSE
    )
  }

  ord <- order(id, repeated)
  y <- as.numeric(y[ord])
  id <- id[ord]
  repeated <- repeated[ord]
  covariates <- covariates[ord, , drop = FALSE]

  current_index <- integer(0)
  previous_index <- integer(0)
  intermittent <- FALSE
  id_split <- split(seq_along(id), id)
  for (idx in id_split) {
    if (length(idx) < 2L) next
    missing_i <- is.na(y[idx])
    if (any(missing_i)) {
      first_missing <- which(missing_i)[1L]
      if (first_missing < length(idx) && any(!missing_i[(first_missing + 1L):length(idx)])) {
        intermittent <- TRUE
      }
    }
    for (k in 2:length(idx)) {
      if (!is.na(y[idx[k - 1L]])) {
        previous_index <- c(previous_index, idx[k - 1L])
        current_index <- c(current_index, idx[k])
      }
    }
  }

  if (!length(current_index)) {
    stop(
      "no transition risk set could be constructed; at least one cluster must have two occasions with the previous response observed",
      call. = FALSE
    )
  }

  list(
    missing = as.integer(is.na(y[current_index])),
    lag_response = y[previous_index],
    id = id[current_index],
    repeated = repeated[current_index],
    occasion = repeated[current_index],
    covariates = covariates[current_index, , drop = FALSE],
    covariate_names = colnames(covariates),
    formula = formula,
    intermittent = intermittent
  )
}


select_mcar_covariates <- function(covariates, occasion, lag_response) {
  occasion_factor <- factor(occasion)
  nuisance <- if (nlevels(occasion_factor) > 1L) {
    stats::model.matrix(~ occasion_factor)
  } else {
    matrix(1, nrow = length(occasion), ncol = 1L)
  }
  base <- nuisance
  rank_base <- qr(base)$rank

  keep <- logical(ncol(covariates))
  if (ncol(covariates)) {
    for (j in seq_len(ncol(covariates))) {
      candidate <- cbind(base, covariates[, j, drop = FALSE])
      rank_candidate <- qr(candidate)$rank
      if (rank_candidate > rank_base) {
        keep[j] <- TRUE
        base <- candidate
        rank_base <- rank_candidate
      }
    }
  }

  lag_candidate <- cbind(base, lag_response)
  if (qr(lag_candidate)$rank <= rank_base) {
    stop(
      "the previous response has no estimable effect after adjustment for occasion and the supplied covariates",
      call. = FALSE
    )
  }

  list(
    covariates = covariates[, keep, drop = FALSE],
    kept_names = colnames(covariates)[keep],
    dropped_names = colnames(covariates)[!keep],
    occasion_factor = occasion_factor
  )
}


fit_mcar_binary_model <- function(formula, analysis_data, orstr, control) {
  mcar_id <- analysis_data[["mcar_id"]]
  mcar_repeated <- analysis_data[["mcar_repeated"]]

  fit <- geewa_binary(
    formula = formula,
    link = "logit",
    data = analysis_data,
    id = mcar_id,
    repeated = mcar_repeated,
    control = control,
    orstr = orstr,
    method = "gee"
  )
  if (!isTRUE(fit$converged)) {
    stop("the binary GEE missingness model did not converge", call. = FALSE)
  }
  fit
}


mcar_nested_geer_test <- function(object0,
                                  object1,
                                  test,
                                  cov_type,
                                  pmethod) {
  value <- switch(
    test,
    wald = wald_test(object0, object1, cov_type),
    score = score_test(object0, object1, cov_type),
    `working-wald` = working_wald_test(
      object0, object1, cov_type, pmethod
    ),
    `working-score` = working_score_test(
      object0, object1, cov_type, pmethod
    ),
    `working-lrt` = working_lrt_test(
      object0, object1, cov_type, pmethod
    )
  )
  list(
    statistic = value$test_stat,
    df = as.numeric(value$test_df),
    p_value = value$test_p
  )
}


mcar_zero_test <- function() {
  list(statistic = 0, df = 0, p_value = 1)
}


mcar_coefficient_table <- function(fit, covariance, covariate_map) {
  beta <- stats::coef(fit)
  standard_error <- sqrt(diag(covariance))
  z_value <- as.numeric(beta) / standard_error
  display_names <- names(beta)
  if (length(covariate_map)) {
    for (safe in names(covariate_map)) {
      display_names[display_names == safe] <- covariate_map[[safe]]
    }
  }
  display_names[display_names == "mcar_lag_response"] <- "lag_response"
  display_names <- sub("^mcar_occasion", "occasion", display_names)

  data.frame(
    term = display_names,
    estimate = as.numeric(beta),
    std.error = standard_error,
    statistic = z_value,
    p.value = 2 * stats::pnorm(abs(z_value), lower.tail = FALSE),
    row.names = NULL,
    check.names = FALSE
  )
}


#' Ridout-Style Regression Diagnostic for MCAR in Longitudinal Data
#'
#' @description
#' Assesses whether longitudinal response missingness is associated with the
#' immediately previous observed response and with observed covariates. The
#' binary missingness model is fitted with \code{\link{geewa_binary}} using a
#' logit link. Occasion effects are included as nuisance terms. The same five
#' hypothesis-testing procedures available in \code{anova.geer} can be used to
#' assess response-history dependence, covariate dependence, and their joint
#' contribution.
#'
#' @param object a fitted object of class \code{"geer"}.
#' @param formula optional one-sided formula specifying fully observed
#'   covariates of response missingness, for example
#'   \code{~ treatment + age}. The previous response is added automatically
#'   and should not be included in \code{formula}. If omitted, the right-hand
#'   side of the fitted GEE mean model is used. Terms aliased with the
#'   occasion effects, such as a main linear time effect, are omitted
#'   automatically and reported in \code{dropped_covariates}.
#' @param data optional original data used to fit \code{object}. This is only
#'   needed when the original data cannot be recovered from the fitted object.
#' @param orstr working odds-ratio structure for the binary missingness GEE.
#'   One of \code{"independence"}, \code{"exchangeable"}, or
#'   \code{"unstructured"}. Defaults to \code{"independence"}.
#' @param test hypothesis-testing procedure. One of \code{"wald"},
#'   \code{"score"}, \code{"working-wald"}, \code{"working-score"}, or
#'   \code{"working-lrt"}. These are the same procedures available in
#'   \code{\link{anova.geer}}. Defaults to \code{"wald"}.
#' @param cov_type covariance estimator used for inference and for the
#'   coefficient table. One of \code{"bias-corrected"}, \code{"robust"},
#'   \code{"df-adjusted"}, or \code{"naive"}. Defaults to
#'   \code{"bias-corrected"}.
#' @param pmethod approximation used to compute p-values for the modified
#'   working tests. One of \code{"rao-scott"} or
#'   \code{"satterthwaite"}. It is ignored for \code{"wald"} and
#'   \code{"score"}. Defaults to \code{"rao-scott"}.
#' @param control a \code{\link{geer_control}} object controlling the
#'   binary GEE fit.
#'
#' @details
#' For each subject, the diagnostic considers transitions from occasion
#' \eqn{t-1} to occasion \eqn{t}. A transition enters the risk set only when
#' \eqn{Y_{i,t-1}} is observed. The binary response is
#' \deqn{R_{it}=I(Y_{it}\text{ is missing}),}
#' and the fitted mean model has the form
#' \deqn{\mathrm{logit}\{P(R_{it}=1)\}
#' = \alpha_t + X_{it}^T\gamma + \delta Y_{i,t-1}.}
#' The occasion-specific effects \eqn{\alpha_t} are nuisance parameters.
#' The model is fitted with \code{geewa_binary(..., link = "logit",
#' method = "gee")} using the original subject identifiers.
#'
#' Three tests using the procedure selected by \code{test} are returned in
#' \code{tests}:
#' \itemize{
#'   \item \code{response_history}: tests \eqn{H_0:\delta=0}. Rejection
#'   indicates that missingness depends on the previously observed response and
#'   is evidence against covariate-dependent MCAR/random dropout.
#'   \item \code{covariates}: tests \eqn{H_0:\gamma=0}. Covariate
#'   dependence alone is evidence against strict MCAR relative to the included
#'   covariates but can remain compatible with covariate-dependent MCAR.
#'   \item \code{overall}: tests \eqn{H_0:\gamma=0,\delta=0}.
#' }
#' The primary \code{"htest"} statistic is the \code{response_history}
#' test because this directly assesses covariate-dependent MCAR/random dropout,
#' the condition most relevant to ordinary GEE consistency. The selected
#' procedure is applied consistently to all three hypotheses. The modified
#' working LRT requires \code{orstr = "independence"}, matching the restriction
#' in \code{anova.geer}.
#'
#' The approach is closely related to Ridout's logistic-regression formulation
#' for studying random dropout. Using GEE for the binary indicators allows
#' correlation among repeated missingness indicators within a subject. The
#' procedure is diagnostic: failure to reject does not prove MCAR and does not
#' rule out dependence on unobserved responses (MNAR).
#'
#' With monotone dropout, the risk-set construction contributes each subject up
#' to the first missing response. With intermittent missingness, transitions are
#' included whenever the immediately previous response is observed. A warning is
#' issued in that case because the result should be interpreted as a local
#' missingness-transition diagnostic rather than a pure dropout test.
#'
#' Covariates in \code{formula} must be fully observed on the reconstructed
#' longitudinal data. Completely absent measurement rows cannot be detected
#' unless they are explicitly represented in the supplied data.
#'
#' @return An object of class \code{"htest"}. In addition to the usual
#' components, it contains:
#' \itemize{
#'   \item \code{tests}: data frame containing the response-history,
#'   covariate, and overall tests using the selected procedure.
#'   \item \code{coefficients}: coefficient table for the binary GEE
#'   missingness model.
#'   \item \code{model}: fitted \code{geer} object returned by
#'   \code{geewa_binary}, or \code{NULL} when no missing response occurs in
#'   the risk set.
#'   \item \code{formula}: one-sided covariate formula.
#'   \item \code{dropped_covariates}: model-matrix columns omitted because
#'   they are aliased with occasion effects or earlier covariate columns.
#'   \item \code{missing}, \code{observed}, and \code{transitions}: numbers
#'   of missing outcomes, observed outcomes, and total transitions in the risk
#'   set.
#'   \item \code{clusters}: number of subjects represented in the risk set.
#'   \item \code{intermittent}: whether an observed response occurs after a
#'   missing response for at least one subject.
#'   \item \code{orstr}, \code{test}, \code{cov_type}, and
#'   \code{pmethod}: association, testing, covariance, and working-test
#'   approximation choices used for inference.
#' }
#'
#' @references
#' Rubin, D.B. (1976). Inference and missing data. \emph{Biometrika},
#' \bold{63}, 581--592. \doi{10.1093/biomet/63.3.581}
#'
#' Ridout, M.S. (1991). Testing for random dropouts in repeated measurement
#' data. \emph{Biometrics}, \bold{47}, 1617--1619.
#' \doi{10.2307/2532413}
#'
#' Fitzmaurice, G.M., Heath, A.F. and Clifford, P. (1996). Logistic regression
#' models for binary panel data with attrition. \emph{Journal of the Royal
#' Statistical Society: Series A}, \bold{159}, 249--263.
#' \doi{10.2307/2983172}
#'
#' @seealso \code{\link{little_mcar_test}}, \code{\link{geewa_binary}},
#'   \code{\link{runs_test}}.
#'
#' @examples
#' set.seed(1)
#' id <- rep(seq_len(60), each = 4)
#' time <- rep(seq_len(4), times = 60)
#' treatment <- rep(rep(c(0, 1), each = 30), each = 4)
#' y_complete <- 2 + 0.4 * treatment + 0.2 * time + rnorm(length(id))
#' y <- y_complete
#' for (i in seq_len(60)) {
#'   idx <- which(id == i)
#'   for (j in 2:4) {
#'     p <- plogis(-2.8 + 0.5 * treatment[idx[j]] +
#'       0.45 * y_complete[idx[j - 1]])
#'     if (rbinom(1, 1, p) == 1) {
#'       y[idx[j:4]] <- NA_real_
#'       break
#'     }
#'   }
#' }
#' dat <- data.frame(id, time, treatment, y)
#'
#' fit <- geewa(
#'   y ~ treatment + time,
#'   data = dat,
#'   id = id,
#'   repeated = time,
#'   family = gaussian(),
#'   corstr = "independence"
#' )
#' out <- mcar_logistic_test(fit, formula = ~ treatment, test = "wald")
#' out
#' out$tests
#'
#' ## Generalized score version of the same diagnostic
#' mcar_logistic_test(fit, formula = ~ treatment, test = "score")
#'
#' @export
mcar_logistic_test <- function(object,
                               formula = NULL,
                               data = NULL,
                               orstr = c(
                                 "independence", "exchangeable", "unstructured"
                               ),
                               test = c(
                                 "wald", "score", "working-wald",
                                 "working-score", "working-lrt"
                               ),
                               cov_type = c(
                                 "bias-corrected", "robust", "df-adjusted", "naive"
                               ),
                               pmethod = c("rao-scott", "satterthwaite"),
                               control = geer_control()) {
  object <- check_geer_object(object)
  orstr <- match.arg(orstr)
  opts <- normalize_geer_test_options(
    test = test[1L],
    cov_type = cov_type[1L],
    pmethod = pmethod[1L]
  )
  test <- opts$test
  cov_type <- opts$cov_type
  pmethod <- opts$pmethod
  if (identical(test, "working-lrt") && !identical(orstr, "independence")) {
    stop(
      "the modified working LRT can only be applied to an independence working model",
      call. = FALSE
    )
  }
  covariate_formula <- make_mcar_covariate_formula(object, formula)
  prepared <- reconstruct_mcar_transition_frame(
    object,
    formula = covariate_formula,
    data = data
  )

  if (prepared$intermittent) {
    warning(
      paste0(
        "intermittent response missingness detected; the diagnostic uses only ",
        "transitions for which the immediately previous response is observed ",
        "and should be interpreted as a local missingness-transition diagnostic ",
        "rather than a pure dropout test"
      ),
      call. = FALSE
    )
  }

  missing_no <- as.integer(sum(prepared$missing == 1L))
  observed_no <- as.integer(sum(prepared$missing == 0L))
  transition_no <- as.integer(length(prepared$missing))
  cluster_no <- as.integer(length(unique(prepared$id)))

  empty_tests <- data.frame(
    test = c("response_history", "covariates", "overall"),
    procedure = rep(test, 3L),
    statistic = c(0, 0, 0),
    df = c(0, 0, 0),
    p.value = c(1, 1, 1),
    row.names = NULL,
    check.names = FALSE
  )

  if (missing_no == 0L) {
    return(structure(
      list(
        statistic = c("X-squared" = 0),
        parameter = c(df = 0L),
        p.value = 1,
        method = paste0(
          "Ridout-style response-history diagnostic for MCAR using binary GEE: ",
          format_test_label(test), " test"
        ),
        data.name = "longitudinal response-missingness transitions from fitted geer object",
        alternative = "missingness depends on the previous observed response after adjustment for covariates and occasion",
        tests = empty_tests,
        coefficients = NULL,
        model = NULL,
        formula = covariate_formula,
        dropped_covariates = character(0),
        missing = missing_no,
        observed = observed_no,
        transitions = transition_no,
        clusters = cluster_no,
        intermittent = prepared$intermittent,
        orstr = orstr,
        test = test,
        cov_type = cov_type,
        pmethod = pmethod
      ),
      class = "htest"
    ))
  }
  if (observed_no == 0L) {
    stop("all risk-set outcomes are missing; the missingness model cannot be fitted", call. = FALSE)
  }
  if (cluster_no < 2L) {
    stop("the missingness model requires at least two independent clusters", call. = FALSE)
  }

  selected <- select_mcar_covariates(
    prepared$covariates,
    prepared$occasion,
    prepared$lag_response
  )

  predictor_matrix <- selected$covariates
  safe_names <- paste0("mcar_x", seq_len(ncol(predictor_matrix)))
  analysis_data <- data.frame(
    mcar_missing = prepared$missing,
    mcar_id = prepared$id,
    mcar_repeated = prepared$repeated,
    mcar_occasion = selected$occasion_factor,
    predictor_matrix,
    mcar_lag_response = prepared$lag_response,
    check.names = FALSE
  )
  if (length(safe_names)) {
    start <- 5L
    finish <- start + length(safe_names) - 1L
    names(analysis_data)[start:finish] <- safe_names
  }

  fit_terms <- c(
    if (nlevels(selected$occasion_factor) > 1L) "mcar_occasion" else character(0),
    safe_names,
    "mcar_lag_response"
  )
  fit_formula <- stats::reformulate(fit_terms, response = "mcar_missing")

  missing_fit <- fit_mcar_binary_model(
    fit_formula, analysis_data, orstr, control
  )

  nuisance_terms <- if (nlevels(selected$occasion_factor) > 1L) {
    "mcar_occasion"
  } else {
    character(0)
  }
  response_null_formula <- stats::reformulate(
    c(nuisance_terms, safe_names),
    response = "mcar_missing"
  )
  overall_null_formula <- stats::reformulate(
    nuisance_terms,
    response = "mcar_missing"
  )
  response_null_fit <- fit_mcar_binary_model(
    response_null_formula, analysis_data, orstr, control
  )
  overall_null_fit <- fit_mcar_binary_model(
    overall_null_formula, analysis_data, orstr, control
  )

  response_test <- mcar_nested_geer_test(
    response_null_fit,
    missing_fit,
    test = test,
    cov_type = cov_type,
    pmethod = pmethod
  )
  overall_test <- mcar_nested_geer_test(
    overall_null_fit,
    missing_fit,
    test = test,
    cov_type = cov_type,
    pmethod = pmethod
  )

  if (length(safe_names)) {
    covariate_null_formula <- stats::reformulate(
      c(nuisance_terms, "mcar_lag_response"),
      response = "mcar_missing"
    )
    covariate_null_fit <- fit_mcar_binary_model(
      covariate_null_formula, analysis_data, orstr, control
    )
    covariate_test <- mcar_nested_geer_test(
      covariate_null_fit,
      missing_fit,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod
    )
  } else {
    covariate_null_fit <- NULL
    covariate_test <- mcar_zero_test()
  }

  tests <- data.frame(
    test = c("response_history", "covariates", "overall"),
    procedure = rep(test, 3L),
    statistic = c(
      response_test$statistic,
      covariate_test$statistic,
      overall_test$statistic
    ),
    df = c(response_test$df, covariate_test$df, overall_test$df),
    p.value = c(
      response_test$p_value,
      covariate_test$p_value,
      overall_test$p_value
    ),
    row.names = NULL,
    check.names = FALSE
  )

  beta <- stats::coef(missing_fit)
  covariance <- stats::vcov(missing_fit, cov_type = cov_type)

  covariate_map <- stats::setNames(selected$kept_names, safe_names)
  coefficient_table <- mcar_coefficient_table(
    missing_fit,
    covariance,
    covariate_map
  )

  structure(
    list(
      statistic = c("X-squared" = response_test$statistic),
      parameter = c(df = response_test$df),
      p.value = response_test$p_value,
      method = paste0(
        "Ridout-style response-history diagnostic for MCAR using binary GEE: ",
        format_test_label(test), " test"
      ),
      data.name = "longitudinal response-missingness transitions from fitted geer object",
      alternative = "missingness depends on the previous observed response after adjustment for covariates and occasion",
      tests = tests,
      coefficients = coefficient_table,
      model = missing_fit,
      formula = covariate_formula,
      dropped_covariates = selected$dropped_names,
      missing = missing_no,
      observed = observed_no,
      transitions = transition_no,
      clusters = cluster_no,
      intermittent = prepared$intermittent,
      orstr = orstr,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod
    ),
    class = "htest"
  )
}
