#' @title
#' Jamshidian-Jalal Screening Diagnostic for MCAR
#'
#' @description
#' Performs the covariance-homogeneity diagnostics proposed by Jamshidian and
#' Jalal (2010) for incomplete multivariate data, following the implementation
#' described by Jamshidian, Jalal and Jansen (2014). The procedure groups cases
#' by their original missingness patterns, imputes the incomplete values, and
#' then assesses whether the covariance structure is homogeneous across the
#' pattern groups. It is a screening diagnostic for MCAR and does not condition
#' on the regression model fitted by \code{geer}.
#'
#' @param object a fitted \code{geer} object, a numeric matrix, or a numeric data
#'   frame containing missing values. For a fitted \code{geer} object, the
#'   original response measurements are reconstructed in subject-by-occasion
#'   form, with rows ordered by the sorted cluster identifiers, before rows
#'   omitted by the model fit are removed.
#' @param data optional original data used to fit \code{object}. This is only
#'   used when \code{object} is a \code{geer} fit and is useful when the
#'   original data cannot be recovered from the fitted object.
#' @param method character string selecting the diagnostic. \code{"hawkins"}
#'   uses the modified Hawkins normal-theory test, \code{"nonparametric"} uses
#'   the k-sample Anderson-Darling test, and \code{"auto"} follows the
#'   diagnostic logic of Jamshidian and Jalal: the Hawkins test is examined
#'   first and the nonparametric test distinguishes nonnormality from covariance
#'   heterogeneity when Hawkins rejects. Defaults to \code{"auto"}.
#' @param imputation imputation method used before the diagnostics are
#'   calculated. \code{"distribution-free"} adds resampled complete-case
#'   residuals to best linear predictors of the missing values, in the spirit
#'   of Srivastava and Dolatabadi (2009), with the location vector and
#'   covariance matrix estimated from the complete cases. Following the
#'   \pkg{MissMech} implementation, it requires at least 10 complete cases and
#'   at least \code{2 * p} complete cases among the rows retained after
#'   applying \code{min_pattern_size}, where \code{p} is the number of
#'   variables; otherwise the function warns and uses \code{"normal"}
#'   imputation. The complete-case pattern is therefore unavailable when it has
#'   fewer than \code{min_pattern_size} cases. \code{"normal"} draws from the conditional multivariate normal
#'   distribution using maximum likelihood estimates obtained by EM. Because
#'   normal-theory imputation can inflate the size of the nonparametric test
#'   when the data are not multivariate normal, a warning is issued when
#'   \code{"normal"} is requested with \code{method} other than
#'   \code{"hawkins"}. This warning is not issued when normal-theory
#'   imputation replaces distribution-free imputation because too few complete
#'   cases are available, although the same caution applies. Must not be
#'   supplied together with \code{imputed_data}. Defaults to
#'   \code{"distribution-free"}.
#' @param n_imputations positive integer giving the number of imputed data
#'   sets. The primary statistic, p-value, detailed components and
#'   interpretation are based on the first imputation, as in \pkg{MissMech};
#'   the results for all imputations are returned in \code{imputations}. The
#'   estimates used for imputation are held fixed across imputations. Must be
#'   \code{1} when \code{imputed_data} is supplied. Defaults to \code{1}.
#' @param imputed_data optional numeric matrix or numeric data frame containing
#'   a completed version of the incomplete data, for example obtained by
#'   another imputation method. It must have the same dimensions and row order
#'   as the matrix or data frame supplied in \code{object}, including any rows
#'   with no observed values, contain no missing or non-finite values and
#'   reproduce the observed values. When supplied, no imputation is performed.
#'   This argument is intended for matrix or data-frame input: for a fitted
#'   \code{geer} object the reconstructed subject-by-occasion matrix is not
#'   returned, so a conformable completed data set cannot be guaranteed.
#'   Defaults to \code{NULL}.
#' @param min_pattern_size integer greater than or equal to 2 specifying the
#'   minimum number of cases required for a missingness pattern to be retained.
#'   Defaults to \code{7}, corresponding to omitting patterns with six or fewer
#'   cases, as in the simulations of Jamshidian and Jalal (2010) and the
#'   \pkg{MissMech} default.
#' @param nrep positive integer giving the number of simulated uniform samples
#'   used to approximate the null distribution of a pattern-specific Neyman
#'   smooth statistic when the pattern contains fewer than \code{n_min} cases.
#'   Not used when \code{method = "nonparametric"}. Defaults to \code{10000}.
#' @param n_min integer greater than or equal to 2 specifying the pattern size
#'   from which the chi-squared approximation with four degrees of freedom is
#'   used for the Neyman smooth statistic instead of its simulated null
#'   distribution. Not used when \code{method = "nonparametric"}. Defaults to
#'   \code{30}.
#' @param alpha a single number strictly between 0 and 1 specifying the
#'   significance level used by the automatic diagnostic rule and its
#'   interpretation. With \code{method = "auto"}, the nonparametric diagnostic
#'   is selected when the modified Hawkins test has p-value less than or equal
#'   to \code{alpha}. Defaults to \code{0.05}.
#' @param seed \code{NULL} or a single whole number not exceeding
#'   \code{.Machine$integer.max} in absolute value, specifying the
#'   random-number seed used for imputation and simulated Neyman p-values. The
#'   previous R random-number state is restored on exit. Use \code{NULL} to use
#'   the current random-number stream. Defaults to \code{110}.
#' @param maxit positive integer giving the maximum number of EM iterations for
#'   normal-theory imputation. Defaults to \code{1000}.
#' @param tol a single positive finite number specifying the relative
#'   convergence tolerance for the EM algorithm. Defaults to \code{1e-8}.
#'
#' @details
#' The diagnostic requires missing values. After applying
#' \code{min_pattern_size}, at least two missingness patterns must remain and at
#' least one retained pattern must contain missing values.
#'
#' The implementation follows \pkg{MissMech} (Jamshidian, Jalal and Jansen,
#' 2014), which differs from Jamshidian and Jalal (2010) in three respects.
#' A single imputation method, distribution-free by default, completes the
#' data for both tests, whereas Jamshidian and Jalal (2010) pair normal-theory
#' imputation with the Hawkins test and distribution-free imputation with the
#' nonparametric test. When there are too few complete cases for the
#' distribution-free imputation, normal-theory imputation replaces it, whereas
#' Jamshidian and Jalal (2010, Section 3.2) substitute the maximum likelihood
#' estimates of the location vector and covariance matrix for their
#' complete-case counterparts within the distribution-free procedure, which is
#' what their simulations use. The chi-squared
#' reference is used for the Neyman statistic in patterns with at least
#' \code{n_min} cases, whereas Jamshidian and Jalal (2010) simulate the null
#' distribution for all pattern sizes; setting \code{n_min} larger than the
#' largest retained pattern reproduces their choice.
#'
#' The modified Hawkins component, based on Hawkins (1981), transforms
#' within-pattern Mahalanobis distances to variables that should be Uniform(0, 1) under multivariate
#' normality and homogeneous covariance matrices. A fourth-order Neyman smooth
#' test is applied within each pattern and the pattern-specific p-values are
#' combined by Fisher's method. For patterns with fewer than \code{n_min}
#' cases, the null distribution of the Neyman statistic is simulated once per
#' pattern and reused across imputations, and the p-value is
#' \eqn{(1 + b) / (nrep + 1)}, where \eqn{b} is the number of simulated
#' statistics at least as large as the observed statistic. \pkg{MissMech}
#' instead uses the proportion of simulated statistics exceeding the observed
#' statistic, replacing a zero proportion by \eqn{1 / nrep}.
#'
#' The nonparametric component compares the distributions of the Hawkins
#' transformed distances across missingness-pattern groups using the k-sample
#' Anderson-Darling statistic of Scholz and Stephens (1987). Its p-value is
#' obtained from the standardized statistic with the simulated reference
#' quantiles and smoothing-spline interpolation of the \pkg{kSamples} package
#' (Scholz and Zhu, 2025), which are tabulated for upper-tail probabilities
#' from \eqn{10^{-5}} to \eqn{1 - 10^{-5}}. Outside this range the log-odds
#' are extrapolated linearly and the \code{extrapolated} component of
#' \code{nonparametric} is \code{TRUE}. When the standardized statistic lies
#' above the tabulated quantiles, the p-value is below \eqn{10^{-5}} and
#' indicates strong evidence against homogeneity, but its magnitude is
#' unreliable; when it lies below them, the p-value exceeds
#' \eqn{1 - 10^{-5}}. \pkg{MissMech} interpolates between the
#' five critical values tabulated by Scholz and Stephens (1987), so p-values
#' can differ from \pkg{MissMech}, particularly below 0.01. This component is
#' intended to reduce sensitivity to departures from multivariate normality.
#'
#' With \code{method = "auto"}, both components are calculated. If the Hawkins
#' test does not reject, there is no evidence against its joint normality and
#' homoscedasticity null. If Hawkins rejects but the nonparametric test does not,
#' the result is consistent with nonnormality rather than covariance
#' heterogeneity. Rejection by the nonparametric test provides evidence against
#' covariance homogeneity and therefore against MCAR in the Jamshidian-Jalal
#' framework.
#'
#' The procedure assumes that, apart from the missingness-pattern grouping, the
#' cases come from a common population. Known substantive groups with genuinely
#' different covariance matrices can therefore cause rejection even when
#' missingness is MCAR. Conversely, supplying complete data through
#' \code{imputed_data}, together with an incomplete copy in which each known
#' group is assigned its own artificial missingness pattern, gives a test of
#' covariance homogeneity across those groups (Jamshidian, Jalal and Jansen,
#' 2014, Example 5).
#'
#' Both components treat the case-level Hawkins statistics within a pattern as
#' approximately independent. Jamshidian and Jalal (2010, Section 5) found this
#' adequate for patterns with at least four cases, and Jamshidian, Jalal and
#' Jansen (2014, Section 3.1) report that retaining patterns with four or more
#' cases works well. Smaller values of \code{min_pattern_size} retain more
#' cases at the cost of less reliable reference distributions.
#'
#' With \code{n_imputations > 1}, the tests are repeated on each completed data
#' set. Jamshidian and Jalal (2010, Section 6) recommend examining the
#' variability of the resulting p-values, for example the proportion of
#' imputations in which the test rejects, and the pattern-specific Neyman
#' p-values and Anderson-Darling contributions, to identify patterns whose
#' covariance structure differs from the rest. These results are exploratory
#' and are not pooled into a single test, because the p-values obtained from
#' different imputations are not independent. Distribution-free imputation
#' resamples the complete-case residuals afresh for each imputed data set while
#' holding the complete-case location vector and covariance matrix fixed, which
#' is the multiple-imputation method Jamshidian and Jalal (2010, Section 6)
#' adopt. Normal-theory imputation likewise holds the maximum likelihood
#' estimates fixed and redraws only the conditional normal deviates, as in
#' \pkg{MissMech}; Jamshidian and Jalal (2010, Section 6) also describe drawing
#' those estimates from their asymptotic distribution at each imputation, which
#' would additionally reflect estimation uncertainty and is not implemented
#' here.
#'
#' This procedure ignores the marginal regression structure and should be used
#' as a screening diagnostic. In longitudinal \code{geer} analyses,
#' \code{\link{mcar_logistic_test}} provides the complementary model-aware
#' diagnostic of whether observed covariates or previous responses predict
#' missingness. Failure to reject either component does not prove MCAR.
#'
#' @return An object inheriting from \code{"htest"} with the primary statistic
#'   and p-value for the selected diagnostic, plus the following additional
#'   components:
#'   \describe{
#'     \item{tests}{A data frame containing the Hawkins and/or nonparametric
#'       statistics and p-values for the first imputation.}
#'     \item{hawkins}{Detailed results from the modified Hawkins test for the
#'       first imputation. With \code{method = "nonparametric"}, only the
#'       transformed statistics \code{f.values} and \code{uniform.values},
#'       which the nonparametric test uses, the pooled covariance matrix and
#'       the pattern means are calculated; the statistic, degrees of freedom
#'       and p-values are \code{NA}.}
#'     \item{nonparametric}{Detailed results from the k-sample
#'       Anderson-Darling test for the first imputation, when calculated. Its
#'       \code{group.statistics} are the pattern contributions
#'       \eqn{T_i / N} to the statistic \eqn{T} of Scholz and Stephens (1987),
#'       so that they sum to \eqn{T}; these are the values reported as
#'       \eqn{T_i} by \pkg{MissMech}. Its \code{extrapolated} component
#'       indicates whether the p-value was extrapolated beyond the reference
#'       quantiles.}
#'     \item{imputations}{A list with components \code{tests}, a data frame
#'       with one row per imputation containing the Hawkins and nonparametric
#'       statistics and p-values; \code{hawkins.group.p.values}, a matrix of
#'       pattern-specific Neyman p-values with one row per imputation; and
#'       \code{nonparametric.group.statistics}, a matrix of pattern
#'       contributions to the Anderson-Darling statistic with one row per
#'       imputation. Quantities that were not calculated are \code{NA} or
#'       \code{NULL}.}
#'     \item{n.imputations}{Number of completed data sets analyzed.}
#'     \item{patterns}{The retained missingness patterns, with 1 denoting a
#'       missing value and 0 an observed value.}
#'     \item{pattern.counts}{Numbers of cases in the retained patterns.}
#'     \item{omitted.patterns}{Patterns omitted because they did not meet
#'       \code{min_pattern_size}.}
#'     \item{n}{Number of rows retained for the diagnostic.}
#'     \item{p}{Number of variables or repeated measurements tested.}
#'     \item{imputation.requested,imputation.used}{Requested and actual
#'       imputation methods; both are \code{"user-supplied"} when
#'       \code{imputed_data} is supplied.}
#'     \item{complete.cases}{Number of complete cases among the rows retained
#'       after applying \code{min_pattern_size}.}
#'     \item{imputed.data}{The first completed data set used by the tests.}
#'     \item{location}{Estimated location vector used for imputation, or
#'       \code{NULL} when \code{imputed_data} is supplied.}
#'     \item{covariance}{Estimated covariance matrix used for imputation, or
#'       \code{NULL} when \code{imputed_data} is supplied.}
#'     \item{em.iterations}{Number of EM iterations used for normal-theory
#'       imputation; \code{NA} for distribution-free or user-supplied
#'       imputation.}
#'     \item{em.converged}{Whether the EM algorithm converged for normal-theory
#'       imputation; \code{NA} for distribution-free or user-supplied
#'       imputation.}
#'     \item{method.requested}{Diagnostic requested through \code{method}.}
#'     \item{selected.test}{Diagnostic supplying the primary statistic and
#'       p-value.}
#'     \item{alpha}{Significance level used by the automatic diagnostic rule.}
#'     \item{interpretation}{A concise interpretation of the first imputation
#'       using \code{alpha}.}
#'   }
#'
#' @references
#' Hawkins, D.M. (1981) A new test for multivariate normality and
#' homoscedasticity. \emph{Technometrics}, \bold{23}, 105--110.
#'
#' Jamshidian, M. and Jalal, S. (2010) Tests of homoscedasticity, normality,
#' and missing completely at random for incomplete multivariate data.
#' \emph{Psychometrika}, \bold{75}, 649--674.
#'
#' Jamshidian, M., Jalal, S. and Jansen, C. (2014) MissMech: An R package for
#' testing homoscedasticity, multivariate normality, and missing completely at
#' random (MCAR). \emph{Journal of Statistical Software}, \bold{56}, 1--31.
#'
#' Scholz, F.W. and Stephens, M.A. (1987) K-sample Anderson-Darling tests.
#' \emph{Journal of the American Statistical Association}, \bold{82}, 918--924.
#'
#' Scholz, F. and Zhu, A. (2025) kSamples: K-sample rank tests and their
#' combinations. R package version 1.2-12.
#' \url{https://CRAN.R-project.org/package=kSamples}
#'
#' Srivastava, M.S. and Dolatabadi, M. (2009) Multiple imputation and other
#' resampling schemes for imputing missing observations. \emph{Journal of
#' Multivariate Analysis}, \bold{100}, 1919--1937.
#'
#' @seealso \code{\link{mcar_little_test}}, \code{\link{mcar_logistic_test}}
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(240), ncol = 3)
#' x[41:60, 3] <- NA
#' x[61:80, 2:3] <- NA
#'
#' mcar_homoscedasticity_test(
#'   x,
#'   method = "nonparametric"
#' )
#'
#' out <- mcar_homoscedasticity_test(
#'   x,
#'   method = "nonparametric",
#'   n_imputations = 5
#' )
#' out$imputations$tests
#'
#' @export
mcar_homoscedasticity_test <- function(
    object,
    data = NULL,
    method = c("auto", "nonparametric", "hawkins"),
    imputation = c("distribution-free", "normal"),
    n_imputations = 1L,
    imputed_data = NULL,
    min_pattern_size = 7L,
    nrep = 10000L,
    n_min = 30L,
    alpha = 0.05,
    seed = 110L,
    maxit = 1000L,
    tol = 1e-8) {
  data_name <- deparse1(substitute(object))
  imputation_supplied <- !missing(imputation)
  method <- match.arg(method)
  imputation <- match.arg(imputation)

  n_imputations <- check_integer_at_least(n_imputations, "n_imputations")
  min_pattern_size <- check_integer_at_least(min_pattern_size, "min_pattern_size", lower = 2L)
  nrep <- check_integer_at_least(nrep, "nrep")
  n_min <- check_integer_at_least(n_min, "n_min", lower = 2L)
  check_probability_open(alpha, "alpha")
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
       abs(seed) > .Machine$integer.max || seed != round(seed))) {
    stop(
      "'seed' must be NULL or a single whole number not exceeding .Machine$integer.max in absolute value",
      call. = FALSE
    )
  }
  maxit <- check_integer_at_least(maxit, "maxit")
  if (!is_positive_scalar(tol)) {
    stop("'tol' must be a single positive finite number", call. = FALSE)
  }

  user_imputed <- !is.null(imputed_data)
  if (user_imputed && imputation_supplied) {
    stop("'imputation' must not be supplied when 'imputed_data' is supplied", call. = FALSE)
  }
  if (user_imputed && n_imputations != 1L) {
    stop("'n_imputations' must be 1 when 'imputed_data' is supplied", call. = FALSE)
  }

  if (inherits(object, "geer")) {
    object <- check_geer_object(object)
    x <- extract_geer_mcar_matrix(object, data = data)
    data_name <- "repeated response from fitted geer object"
  } else {
    if (!is.null(data)) {
      stop("'data' is only used when 'object' is a fitted 'geer' object", call. = FALSE)
    }
    x <- object
  }
  x <- check_mcar_homoscedasticity_matrix(x)
  if (!anyNA(x)) {
    stop(
      "the Jamshidian-Jalal diagnostic requires missing values",
      call. = FALSE
    )
  }
  if (user_imputed) {
    imputed_data <- check_mcar_homoscedasticity_imputed(imputed_data, x)
  }

  pattern_info <- compute_pattern_information(x, min_pattern_size = min_pattern_size)
  x_used <- check_mcar_homoscedasticity_matrix(pattern_info$x)
  attr(x_used, "row.index") <- NULL
  attr(x_used, "original.nrow") <- NULL
  if (!anyNA(x_used)) {
    stop(
      "the Jamshidian-Jalal diagnostic requires missing values in at least one retained pattern",
      call. = FALSE
    )
  }

  if (!user_imputed && identical(imputation, "normal") &&
      !identical(method, "hawkins")) {
    warning(
      paste0(
        "normal-theory imputation can inflate the size of the nonparametric test ",
        "when the data are not multivariate normal; ",
        "'imputation = \"distribution-free\"' is recommended"
      ),
      call. = FALSE
    )
  }

  run_hawkins <- !identical(method, "nonparametric")
  run_nonparametric <- method %in% c("auto", "nonparametric")

  result <- run_with_seed(seed, {
    if (user_imputed) {
      setup <- list(
        mu = NULL,
        sigma = NULL,
        iterations = NA_integer_,
        converged = NA,
        requested = "user-supplied",
        used = "user-supplied",
        n_complete = sum(stats::complete.cases(x_used))
      )
      completed <- list(imputed_data[pattern_info$rows, , drop = FALSE])
    } else {
      setup <- build_imputation_setup(
        x_used,
        imputation = imputation,
        maxit = maxit,
        tol = tol
      )
      completed <- lapply(
        seq_len(n_imputations),
        function(i) draw_imputation(x_used, setup)
      )
    }

    neyman_nulls <- if (run_hawkins) {
      simulate_neyman_nulls(pattern_info$group_counts, nrep = nrep, n_min = n_min)
    } else {
      NULL
    }

    analyses <- lapply(
      completed,
      function(completed_data) {
        hawkins <- compute_hawkins_test(
          completed = completed_data,
          group = pattern_info$group,
          group_counts = pattern_info$group_counts,
          neyman_nulls = neyman_nulls,
          test_uniformity = run_hawkins
        )
        nonparametric <- NULL
        if (run_nonparametric) {
          nonparametric <- compute_nonparametric_test(
            hawkins = hawkins,
            group = pattern_info$group,
            group_counts = pattern_info$group_counts
          )
        }
        list(hawkins = hawkins, nonparametric = nonparametric)
      }
    )

    list(setup = setup, completed = completed, analyses = analyses)
  })

  first <- result$analyses[[1L]]
  hawkins <- first$hawkins
  nonparametric <- first$nonparametric

  if (identical(method, "hawkins")) {
    statistic <- c("Fisher chi-squared" = hawkins$statistic)
    parameter <- c(df = hawkins$parameter)
    p_value <- hawkins$p.value
    method_label <- "Jamshidian-Jalal modified Hawkins MCAR screening diagnostic"
    selected_test <- "hawkins"
  } else if (identical(method, "nonparametric")) {
    statistic <- c("Anderson-Darling" = nonparametric$statistic)
    parameter <- NULL
    p_value <- nonparametric$p.value
    method_label <- "Jamshidian-Jalal nonparametric MCAR screening diagnostic"
    selected_test <- "nonparametric"
  } else if (hawkins$p.value > alpha) {
    statistic <- c("Fisher chi-squared" = hawkins$statistic)
    parameter <- c(df = hawkins$parameter)
    p_value <- hawkins$p.value
    method_label <- "Jamshidian-Jalal automatic MCAR screening diagnostic (modified Hawkins)"
    selected_test <- "hawkins"
  } else {
    statistic <- c("Anderson-Darling" = nonparametric$statistic)
    parameter <- NULL
    p_value <- nonparametric$p.value
    method_label <- "Jamshidian-Jalal automatic MCAR screening diagnostic (nonparametric)"
    selected_test <- "nonparametric"
  }

  tests <- data.frame(
    test = character(),
    statistic = numeric(),
    df = numeric(),
    p.value = numeric(),
    stringsAsFactors = FALSE
  )
  if (run_hawkins) {
    tests <- rbind(
      tests,
      data.frame(
        test = "hawkins",
        statistic = hawkins$statistic,
        df = hawkins$parameter,
        p.value = hawkins$p.value,
        stringsAsFactors = FALSE
      )
    )
  }
  if (run_nonparametric) {
    tests <- rbind(
      tests,
      data.frame(
        test = "nonparametric",
        statistic = nonparametric$statistic,
        df = NA_real_,
        p.value = nonparametric$p.value,
        stringsAsFactors = FALSE
      )
    )
  }

  pattern_names <- rownames(pattern_info$pattern_matrix)
  n_completed <- length(result$analyses)
  extract_value <- function(test, component) {
    vapply(
      result$analyses,
      function(analysis) {
        if (is.null(analysis[[test]])) NA_real_ else analysis[[test]][[component]]
      },
      numeric(1)
    )
  }
  imputations <- list(
    tests = data.frame(
      imputation = seq_len(n_completed),
      hawkins.statistic = extract_value("hawkins", "statistic"),
      hawkins.p.value = extract_value("hawkins", "p.value"),
      nonparametric.statistic = extract_value("nonparametric", "statistic"),
      nonparametric.p.value = extract_value("nonparametric", "p.value"),
      stringsAsFactors = FALSE
    ),
    hawkins.group.p.values = if (run_hawkins) {
      do.call(
        rbind,
        lapply(result$analyses, function(analysis) analysis$hawkins$group.p.values)
      )
    } else {
      NULL
    },
    nonparametric.group.statistics = if (run_nonparametric) {
      do.call(
        rbind,
        lapply(result$analyses, function(analysis) analysis$nonparametric$group.statistics)
      )
    } else {
      NULL
    }
  )
  if (!is.null(imputations$hawkins.group.p.values)) {
    dimnames(imputations$hawkins.group.p.values) <- list(NULL, pattern_names)
  }
  if (!is.null(imputations$nonparametric.group.statistics)) {
    dimnames(imputations$nonparametric.group.statistics) <- list(NULL, pattern_names)
  }

  out <- list(
    statistic = statistic,
    parameter = parameter,
    p.value = p_value,
    method = method_label,
    data.name = data_name,
    alternative = if (identical(selected_test, "hawkins")) {
      "multivariate normality or covariance homogeneity fails"
    } else {
      "covariance structure differs across missingness-pattern groups"
    },
    tests = tests,
    hawkins = hawkins,
    nonparametric = nonparametric,
    imputations = imputations,
    n.imputations = n_completed,
    patterns = pattern_info$pattern_matrix,
    pattern.counts = stats::setNames(pattern_info$group_counts, pattern_names),
    omitted.patterns = pattern_info$omitted_patterns,
    n = nrow(x_used),
    p = ncol(x_used),
    imputation.requested = result$setup$requested,
    imputation.used = result$setup$used,
    complete.cases = result$setup$n_complete,
    imputed.data = result$completed[[1L]],
    location = result$setup$mu,
    covariance = result$setup$sigma,
    em.iterations = result$setup$iterations,
    em.converged = result$setup$converged,
    method.requested = method,
    selected.test = selected_test,
    alpha = alpha,
    interpretation = interpret_homoscedasticity_test(
      method = method,
      hawkins = hawkins,
      nonparametric = nonparametric,
      alpha = alpha
    )
  )
  class(out) <- c("mcar_homoscedasticity_test", "htest")
  out
}
