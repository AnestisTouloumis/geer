#' @title
#' Ridout-Style Regression Diagnostic for MCAR in Longitudinal Data
#'
#' @description
#' Assesses whether longitudinal response missingness is associated with the
#' immediately preceding response, when observed, and with observed covariates. The
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
#'   automatically and reported in \code{dropped_covariates}. Defaults to
#'   \code{NULL}.
#' @param data optional original data used to fit \code{object}. This is only
#'   needed when the original data cannot be recovered from the fitted object.
#'   Defaults to \code{NULL}.
#' @param orstr working odds-ratio structure for the binary missingness GEE.
#'   One of \code{"independence"}, \code{"exchangeable"}, or
#'   \code{"unstructured"}. Defaults to \code{"independence"}.
#' @param test hypothesis-testing procedure. One of \code{"wald"},
#'   \code{"score"}, \code{"working-wald"}, \code{"working-score"}, or
#'   \code{"working-lrt"}. These are the same procedures available in
#'   \code{\link{anova.geer}}. Defaults to \code{"wald"}.
#' @param cov_type covariance estimator used for inference and for the
#'   coefficient table. One of \code{"bias-corrected"}, \code{"robust"},
#'   \code{"df-adjusted"}, \code{"jackknife"}, or \code{"naive"}. Defaults to
#'   \code{"bias-corrected"}.
#' @param pmethod approximation used to compute p-values for the modified
#'   working tests. One of \code{"rao-scott"} or
#'   \code{"satterthwaite"}. It is ignored for \code{"wald"} and
#'   \code{"score"}. Defaults to \code{"rao-scott"}.
#' @param control a list returned by \code{\link{geer_control}} controlling
#'   the binary GEE fit. Defaults to \code{geer_control()}.
#'
#' @details
#' For each cluster, the diagnostic considers transitions from occasion
#' \eqn{t-1} to occasion \eqn{t}. A transition enters the risk set only when
#' \eqn{Y_{i,t-1}} is observed. The binary response is
#' \deqn{R_{it}=I(Y_{it}\text{ is missing}),}
#' and the fitted mean model has the form
#' \deqn{\mathrm{logit}\{P(R_{it}=1)\}
#' = \alpha_t + X_{it}^T\gamma + \delta Y_{i,t-1}.}
#' The occasion-specific effects \eqn{\alpha_t} are nuisance parameters.
#' The model is fitted with \code{geewa_binary(..., link = "logit",
#' method = "gee")} using the original cluster identifiers.
#' If a score-based procedure is combined with \code{cov_type = "jackknife"},
#' the score and model-based information are evaluated under the corresponding
#' null model, while the covariance component is obtained from full
#' leave-one-cluster refits of the larger missingness model.
#'
#' Three tests using the procedure selected by \code{test} are returned in
#' \code{tests}:
#' \itemize{
#'   \item \code{response_history}: tests \eqn{H_0:\delta=0}. Rejection
#'   indicates that missingness depends on the immediately preceding response and
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
#' correlation among repeated missingness indicators within a cluster. The
#' procedure is diagnostic: failure to reject does not prove MCAR and does not
#' rule out dependence on unobserved responses (MNAR).
#'
#' With monotone dropout, the risk-set construction contributes each cluster up
#' to the first missing response. With intermittent missingness, transitions are
#' included whenever the immediately preceding response is observed. A warning is
#' issued in that case because the result should be interpreted as a local
#' missingness-transition diagnostic rather than a pure dropout test.
#'
#' Covariates in \code{formula} must be fully observed on the reconstructed
#' longitudinal data. Completely absent measurement rows cannot be detected
#' unless they are explicitly represented in the supplied data.
#'
#' @return An object of class \code{c("geer_htest", "htest")}. In addition to the usual
#' components, it contains:
#' \itemize{
#'   \item \code{tests}: data frame containing the response-history,
#'   covariate, and overall tests using the selected procedure.
#'   \item \code{coefficients}: coefficient table for the binary GEE
#'   missingness model, or \code{NULL} when no missing response occurs in the
#'   risk set.
#'   \item \code{model}: fitted \code{geer} object returned by
#'   \code{geewa_binary}, or \code{NULL} when no missing response occurs in
#'   the risk set.
#'   \item \code{formula}: one-sided covariate formula.
#'   \item \code{dropped_covariates}: model-matrix columns omitted because
#'   they are aliased with occasion effects or earlier covariate columns.
#'   \item \code{missing}, \code{observed}, and \code{transitions}: numbers
#'   of missing outcomes, observed outcomes, and total transitions in the risk
#'   set.
#'   \item \code{clusters}: number of clusters represented in the risk set.
#'   \item \code{intermittent}: whether an observed response occurs after a
#'   missing response for at least one cluster.
#'   \item \code{orstr}, \code{test}, \code{cov_type}, and
#'   \code{pmethod}: association, testing, covariance, and working-test
#'   approximation choices used for inference.
#' }
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
#' @seealso \code{\link{mcar_little_test}}, \code{\link{geewa_binary}},
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
      "the modified working LRT can only be applied to an independence working model",
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
      class = c("geer_htest", "htest")
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
    covariate_null_formula <- stats::reformulate(
      c(nuisance_terms, "mcar_lag_response"),
      response = "mcar_missing"
    )
    covariate_null_fit <- fit_mcar_binary_model(
      covariate_null_formula, analysis_data, orstr, control
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
    class = c("geer_htest", "htest")
  )
}


