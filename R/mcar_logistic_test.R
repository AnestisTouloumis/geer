#' @title
#' Logistic Regression Test of Missing Completely at Random
#'
#' @description
#' Tests whether the missingness of a longitudinal response depends on the
#' previous response or on observed covariates, using a logistic regression in
#' the spirit of Ridout (1991) and Fitzmaurice, Heath and Clifford (1996) fitted
#' by generalized estimating equations. MCAR is defined as in Rubin (1976).
#'
#' @param object a fitted model object of class \code{"geer"}.
#' @param formula optional one-sided formula of fully observed covariates of
#'   missingness, for example \code{~ treatment + age}. The previous response is
#'   added automatically and must not be included. Terms aliased with the
#'   occasion effects are dropped and reported in \code{dropped_covariates}.
#'   Defaults to \code{NULL}, which uses the right-hand side of the mean model
#'   of \code{object}.
#' @param data optional data frame containing the data used to fit
#'   \code{object}. Used only when \code{object} is a fitted model object and
#'   the data cannot be recovered from it. Defaults to \code{NULL}.
#' @param orstr character string specifying the working odds-ratio structure of
#'   the binary missingness model: \code{"independence"}, \code{"exchangeable"},
#'   or \code{"unstructured"}. Defaults to \code{"independence"}.
#' @param test character string specifying the hypothesis test: Wald
#'   (\code{"wald"}), generalized score (\code{"score"}), modified working Wald
#'   (\code{"working-wald"}), modified working score (\code{"working-score"}),
#'   or modified working likelihood ratio (\code{"working-lrt"}). Defaults to
#'   \code{"wald"}.
#' @param cov_type character string specifying the covariance matrix estimator
#'   used for inference: \code{"bias-corrected"}, \code{"robust"},
#'   \code{"df-adjusted"}, \code{"jackknife"}, or \code{"naive"}. See
#'   \code{\link{vcov.geer}}. Defaults to \code{"bias-corrected"}.
#' @param pmethod character string specifying the approximation to the null
#'   distribution of the modified working tests: Rao-Scott (\code{"rao-scott"})
#'   or Satterthwaite (\code{"satterthwaite"}). Ignored for the Wald and score
#'   tests. Defaults to \code{"rao-scott"}.
#' @param control list of control parameters as returned by
#'   \code{\link{geer_control}} for the binary GEE fit. Defaults to
#'   \code{geer_control()}.
#'
#' @details
#' For each cluster \eqn{i} and occasion \eqn{t}, let
#' \eqn{R_{it} = I(Y_{it} \text{ is missing})}{R_it = I(Y_it is missing)}. The transitions with
#' \eqn{Y_{i,t-1}} observed form the risk set, on which the model
#' \deqn{\mathrm{logit}\{P(R_{it} = 1)\} = \alpha_t + X_{it}^\top\gamma +
#' \delta Y_{i,t-1}}{logit P(R_it = 1) = alpha_t + X_it' gamma + delta Y_i,t-1}
#' is fitted by \code{\link{geewa_binary}} with \code{method = "gee"}, so that
#' the missingness indicators of a cluster may be correlated. The occasion
#' effects \eqn{\alpha_t}{alpha_t} are nuisance parameters. The component \code{tests}
#' contains three tests, all computed by the procedure selected by
#' \code{test}:
#' \describe{
#'   \item{\code{response_history}}{\eqn{H_0: \delta = 0}{H0: delta = 0}. Rejection indicates
#'     that missingness depends on the previous response, which is evidence
#'     against covariate-dependent MCAR. This is the primary test.}
#'   \item{\code{covariates}}{\eqn{H_0: \gamma = 0}{H0: gamma = 0}. Rejection is evidence
#'     against MCAR but is compatible with covariate-dependent MCAR.}
#'   \item{\code{overall}}{\eqn{H_0: \gamma = 0, \delta = 0}{H0: gamma = 0, delta = 0}.}
#' }
#' The modified working likelihood ratio test requires
#' \code{orstr = "independence"}. With score tests and
#' \code{cov_type = "jackknife"}, the covariance component is computed from
#' the larger model.
#'
#' Under monotone dropout each cluster contributes up to its first missing
#' response. With intermittent missingness a warning is issued, because the
#' result is then a local transition diagnostic rather than a dropout test.
#' Covariates in \code{formula} must be fully observed, and measurement rows
#' that are absent from the data cannot be detected: when no response is
#' missing but some clusters have no rows for some occasions, a warning is
#' issued. The test is a
#' diagnostic: failure to reject does not establish MCAR and does not rule
#' out dependence on unobserved responses.
#'
#' @return
#' An object of class \code{c("geer_htest", "htest")}. In addition to the
#' standard components, the object contains:
#' \item{tests}{a data frame containing the response-history, covariate, and
#'   overall tests using the selected procedure.}
#' \item{coefficients}{the coefficient table of the binary GEE missingness
#'   model, or \code{NULL} when no missing response occurs in the risk set.}
#' \item{model}{the fitted model object of class \code{"geer"} returned by
#'   \code{\link{geewa_binary}}, or \code{NULL} when no missing response
#'   occurs in the risk set.}
#' \item{formula}{the one-sided covariate formula.}
#' \item{dropped_covariates}{the model-matrix columns omitted because they are
#'   aliased with the occasion effects or with earlier covariate columns.}
#' \item{missing, observed, transitions}{the numbers of missing outcomes,
#'   observed outcomes, and transitions in the risk set.}
#' \item{clusters}{the number of clusters represented in the risk set.}
#' \item{intermittent}{whether an observed response occurs after a missing
#'   response for at least one cluster.}
#' \item{orstr, test, cov_type, pmethod}{the association structure, testing
#'   procedure, covariance estimator, and working-test approximation used for
#'   inference.}
#'
#' @references
#' Fitzmaurice, G.M., Heath, A.F. and Clifford, P. (1996) Logistic regression
#' models for binary panel data with attrition. \emph{Journal of the Royal
#' Statistical Society: Series A}, \bold{159}, 249--263.
#'
#' Ridout, M.S. (1991) Testing for random dropouts in repeated measurement
#' data. \emph{Biometrics}, \bold{47}, 1617--1619.
#'
#' Rubin, D.B. (1976) Inference and missing data. \emph{Biometrika},
#' \bold{63}, 581--592.
#'
#' @seealso \code{\link{mcar_little_test}},
#'   \code{\link{mcar_homoscedasticity_test}}, \code{\link{geewa}},
#'   \code{\link{geewa_binary}}.
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
#'   formula = y ~ treatment + time,
#'   family = gaussian(),
#'   data = dat,
#'   id = id,
#'   repeated = time,
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
                               orstr = c("independence", "exchangeable",
                                         "unstructured"),
                               test = c("wald", "score", "working-wald",
                                        "working-score", "working-lrt"),
                               cov_type = c("bias-corrected", "robust",
                                            "df-adjusted", "jackknife",
                                            "naive"),
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
      "the modified working LRT requires an independence working structure",
      call. = FALSE
    )
  }
  covariate_formula <- build_mcar_covariate_formula(object, formula)
  prepared <- reconstruct_mcar_transition_frame(
    object,
    formula = covariate_formula,
    data = data
  )

  if (prepared$intermittent) {
    warning(
      paste0(
        "intermittent missingness detected: only transitions with an observed ",
        "previous response are used, so the result is a local transition ",
        "diagnostic rather than a dropout test"
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
    if (isTRUE(prepared$incomplete_clusters)) {
      warning(
        paste0(
          "no response is missing, but some clusters lack rows for some ",
          "occasions; dropout recorded by omitted rows cannot be detected ",
          "unless those rows are added with missing responses"
        ),
        call. = FALSE
      )
    }
    return(structure(
      list(
        statistic = c("X-squared" = 0),
        parameter = c(df = 0L),
        p.value = 1,
        method = paste0(
          "Ridout-style response-history diagnostic for MCAR using binary ",
          "GEE: ",
          format_test_label(test), " test"
        ),
        data.name = paste0(
          "longitudinal response-missingness transitions from fitted geer ",
          "object"
        ),
        alternative = paste0(
          "missingness depends on the previous observed response after ",
          "adjustment for covariates and occasion"
        ),
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
      class = c("geer_htest", "htest")
    ))
  }
  if (observed_no == 0L) {
    stop(
      paste0(
        "all risk-set outcomes are missing; the missingness model cannot be ",
        "fitted"
      ),
      call. = FALSE
    )
  }
  if (cluster_no < 2L) {
    stop(
      "the missingness model requires at least two independent clusters",
      call. = FALSE
    )
  }

  selected <- select_mcar_covariates(
    prepared$covariates,
    prepared$occasion,
    prepared$lag_response
  )

  predictor_matrix <- selected$covariates
  safe_names <- if (ncol(predictor_matrix) > 0L) {
    paste0("mcar_x", seq_len(ncol(predictor_matrix)))
  } else {
    character(0)
  }
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
    if (nlevels(selected$occasion_factor) > 1L) {
      "mcar_occasion"
    } else {
      character(0)
    },
    safe_names,
    "mcar_lag_response"
  )
  fit_formula <- stats::reformulate(fit_terms, response = "mcar_missing")

  missing_fit <- fit_mcar_binary_model(
    fit_formula,
    analysis_data,
    orstr,
    control
  )

  nuisance_terms <- if (nlevels(selected$occasion_factor) > 1L) {
    "mcar_occasion"
  } else {
    character(0)
  }
  response_null_formula <- build_mcar_null_formula(
    c(nuisance_terms, safe_names)
  )
  overall_null_formula <- build_mcar_null_formula(nuisance_terms)
  response_null_fit <- fit_mcar_binary_model(
    response_null_formula,
    analysis_data,
    orstr,
    control
  )
  overall_null_fit <- fit_mcar_binary_model(
    overall_null_formula,
    analysis_data,
    orstr,
    control
  )

  response_test <- compute_mcar_nested_test(
    response_null_fit,
    missing_fit,
    test = test,
    cov_type = cov_type,
    pmethod = pmethod
  )
  overall_test <- compute_mcar_nested_test(
    overall_null_fit,
    missing_fit,
    test = test,
    cov_type = cov_type,
    pmethod = pmethod
  )

  if (length(safe_names)) {
    covariate_null_formula <- build_mcar_null_formula(
      c(nuisance_terms, "mcar_lag_response")
    )
    covariate_null_fit <- fit_mcar_binary_model(
      covariate_null_formula,
      analysis_data,
      orstr,
      control
    )
    covariate_test <- compute_mcar_nested_test(
      covariate_null_fit,
      missing_fit,
      test = test,
      cov_type = cov_type,
      pmethod = pmethod
    )
  } else {
    covariate_null_fit <- NULL
    covariate_test <- build_mcar_zero_test()
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

  covariance <- stats::vcov(missing_fit, cov_type = cov_type)

  covariate_map <- if (length(safe_names)) {
    stats::setNames(selected$kept_names, safe_names)
  } else {
    character(0)
  }
  coefficient_table <- build_mcar_coefficient_table(
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
      data.name = paste0(
        "longitudinal response-missingness transitions from fitted geer object"
      ),
      alternative = paste0(
        "missingness depends on the previous observed response after ",
        "adjustment for covariates and occasion"
      ),
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
    class = c("geer_htest", "htest")
  )
}
