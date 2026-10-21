#' @title
#' Jamshidian-Jalal Test of Missing Completely at Random
#'
#' @description
#' Tests the homogeneity of the covariance matrices across missingness
#' patterns, a consequence of MCAR, by the procedure of Jamshidian and Jalal
#' (2010) as implemented in \pkg{MissMech} (Jamshidian, Jalal and Jansen,
#' 2014).
#'
#' @param object a fitted model object of class \code{"geer"}, or a numeric
#'   matrix or data frame with missing values. For a fitted model object, the
#'   response is arranged with one row per cluster and one column per occasion.
#' @param data optional data frame containing the data used to fit
#'   \code{object}. Used only when \code{object} is a fitted model object and
#'   the data cannot be recovered from it. Defaults to \code{NULL}.
#' @param method character string specifying the diagnostic: the modified
#'   Hawkins test (\code{"hawkins"}), the k-sample Anderson-Darling test
#'   (\code{"nonparametric"}), or the two-stage rule of Jamshidian and Jalal
#'   (2010) (\code{"auto"}). See Details. Defaults to \code{"auto"}.
#' @param imputation character string specifying how the data are completed
#'   before testing: by adding resampled complete-case residuals to the best
#'   linear predictors of the missing values (\code{"distribution-free"}), or by
#'   drawing from the conditional multivariate normal distribution
#'   (\code{"normal"}). See Details. Must not be combined with
#'   \code{imputed_data}. Defaults to \code{"distribution-free"}.
#' @param n_imputations positive integer giving the number of imputed data sets.
#'   The primary results use the first one, as in \pkg{MissMech}, and the
#'   results for all of them are returned in \code{imputations}. Must be 1 when
#'   \code{imputed_data} is supplied. Defaults to \code{1}.
#' @param imputed_data optional completed version of the matrix or data frame
#'   supplied in \code{object}, with the same dimensions and row order, no
#'   missing values, and the observed values unchanged. No imputation is then
#'   performed. Defaults to \code{NULL}.
#' @param min_pattern_size integer, at least 2, giving the minimum number of
#'   cases for a missingness pattern to be retained. Defaults to \code{7}, as in
#'   Jamshidian and Jalal (2010) and \pkg{MissMech}.
#' @param n_replicates positive integer giving the number of simulated samples
#'   for the null distribution of the Neyman statistic in patterns with fewer
#'   than \code{n_min} cases. Defaults to \code{10000}.
#' @param n_min integer, at least 2, giving the pattern size from which the
#'   chi-squared reference with four degrees of freedom is used for the Neyman
#'   statistic. Defaults to \code{30}.
#' @param sig_level a single number strictly between 0 and 1 giving the
#'   significance level of the rule used with \code{method = "auto"}. Defaults
#'   to \code{0.05}.
#' @param seed \code{NULL} or a single integer giving the random number seed
#'   used for imputation and for the simulated Neyman p-values. The previous
#'   random number state is restored on exit, and \code{NULL} uses the current
#'   random number stream. Defaults to \code{110}.
#' @param maxiter positive integer giving the maximum number of iterations of
#'   the EM algorithm used for normal-theory imputation. Defaults to
#'   \code{1000}.
#' @param tolerance a single positive finite number giving the relative
#'   convergence tolerance of the EM algorithm. Defaults to \code{1e-08}.
#'
#' @details
#' The procedure completes the incomplete data by imputation, groups the cases
#' by their original missingness pattern, and tests whether the completed data
#' have a common covariance matrix across patterns. Patterns with fewer than
#' \code{min_pattern_size} cases are omitted; at least two patterns, one of
#' them with missing values, must remain.
#'
#' Two tests are available. The modified Hawkins test (Hawkins, 1981)
#' transforms the within-pattern Mahalanobis distances to variables that are
#' Uniform(0, 1) under multivariate normality and homogeneous covariances,
#' applies a fourth-order Neyman smooth test within each pattern, and combines
#' the pattern p-values by Fisher's method. For patterns with fewer than
#' \code{n_min} cases the null distribution of the Neyman statistic is
#' simulated and the p-value is \eqn{(1 + b)/(n_replicates + 1)}, where
#' \eqn{b} is the number of simulated statistics at least as large as the
#' observed one. The nonparametric test compares the distributions of the
#' transformed distances
#' across patterns with the k-sample Anderson-Darling test of Scholz and
#' Stephens (1987), which is less sensitive to nonnormality. With
#' \code{method = "auto"} both tests are computed: if the Hawkins test does not
#' reject, there is no evidence against MCAR; if it rejects but the
#' nonparametric test does not, the rejection is attributed to nonnormality;
#' rejection by the nonparametric test is evidence against covariance
#' homogeneity and hence against MCAR.
#'
#' Distribution-free imputation, in the spirit of Srivastava and Dolatabadi
#' (2009), requires at least 10, and at least \eqn{2p}, complete cases among the
#' retained rows, where \eqn{p} is the number of variables; otherwise
#' normal-theory imputation is used with a warning. Because normal-theory
#' imputation can inflate the size of the nonparametric test for nonnormal data,
#' requesting it with \code{method} other than \code{"hawkins"} gives a warning.
#' \code{imputed_data} is intended for matrix or data frame input, since the
#' matrix reconstructed from a fitted model object is not returned.
#'
#' The implementation follows \pkg{MissMech} (Jamshidian, Jalal and Jansen,
#' 2014), which differs from Jamshidian and Jalal (2010) in that one imputation
#' method serves both tests, normal-theory imputation replaces
#' distribution-free imputation when there are too few complete cases, and the
#' chi-squared reference is used for the Neyman statistic in patterns with at
#' least \code{n_min} cases. Setting \code{n_min} larger than the largest
#' pattern simulates the reference for all patterns, as in Jamshidian and Jalal
#' (2010). This function differs from \pkg{MissMech} in two respects: the
#' simulated Neyman p-values use \eqn{(1 + b)/(n_replicates + 1)}, and the
#' Anderson-Darling p-value is interpolated from the simulated quantiles of the
#' \pkg{kSamples} package (Scholz and Zhu, 2025) instead of the five critical
#' values of Scholz and Stephens (1987). P-values can therefore differ from
#' \pkg{MissMech}, particularly below 0.01. The quantiles cover upper-tail
#' probabilities from \eqn{10^{-5}} to \eqn{1 - 10^{-5}}; beyond this range the
#' p-value is extrapolated and \code{nonparametric$extrapolated} is
#' \code{TRUE}.
#'
#' With \code{n_imputations > 1} the tests are repeated on each completed data
#' set, holding the estimates used for imputation fixed. Jamshidian and Jalal
#' (2010, Section 6) suggest examining the variability of the p-values and the
#' pattern-specific results to identify patterns whose covariance differs.
#' These results are exploratory and are not pooled.
#'
#' Both tests treat the cases within a pattern as approximately independent,
#' which Jamshidian and Jalal (2010) and Jamshidian, Jalal and Jansen (2014)
#' found adequate for patterns with at least four cases. Known groups with
#' different covariance matrices can cause rejection even under MCAR;
#' conversely, supplying complete data through \code{imputed_data} with each
#' group given its own artificial pattern tests covariance homogeneity across
#' the groups (Jamshidian, Jalal and Jansen, 2014, Example 5).
#'
#' The procedure ignores the marginal regression model and is a screening
#' diagnostic: failure to reject does not establish MCAR. For longitudinal
#' \code{geer} fits, \code{\link{mcar_logistic_test}} assesses whether
#' covariates or previous responses predict missingness.
#'
#' @return
#' An object of class \code{c("geer_htest", "htest")}. The standard
#' components hold the statistic and p-value of the selected test. In
#' addition, the object contains:
#' \item{tests}{a data frame with the statistics and p-values of the tests for
#'   the first imputation.}
#' \item{hawkins}{the detailed results of the Hawkins test for the first
#'   imputation. With \code{method = "nonparametric"}, only the quantities used
#'   by the nonparametric test are computed and the others are \code{NA}.}
#' \item{nonparametric}{the detailed results of the Anderson-Darling test for
#'   the first imputation, if computed, including the pattern contributions
#'   \code{group.statistics}, reported as \eqn{T_i} by \pkg{MissMech}, and the
#'   indicator \code{extrapolated}.}
#' \item{imputations}{a list with the statistics and p-values
#'   (\code{tests}), the pattern Neyman p-values
#'   (\code{hawkins.group.p.values}), and the pattern Anderson-Darling
#'   contributions (\code{nonparametric.group.statistics}) of every
#'   imputation.}
#' \item{n.imputations}{the number of imputed data sets.}
#' \item{patterns, pattern.counts}{the retained missingness patterns, with 1
#'   for a missing value, and their numbers of cases.}
#' \item{omitted.patterns}{the patterns omitted by \code{min_pattern_size}.}
#' \item{n, p}{the numbers of retained rows and of variables.}
#' \item{imputation.requested, imputation.used}{the requested and the actual
#'   imputation method, \code{"user-supplied"} with \code{imputed_data}.}
#' \item{complete.cases}{the number of retained complete cases.}
#' \item{imputed.data}{the first completed data set.}
#' \item{location, covariance}{the estimates used for imputation, \code{NULL}
#'   with \code{imputed_data}.}
#' \item{em.iterations, em.converged}{the number of EM iterations and whether
#'   the EM algorithm converged, \code{NA} unless normal-theory imputation was
#'   used.}
#' \item{method.requested, selected.test}{the requested diagnostic and the
#'   test giving the primary result.}
#' \item{sig.level}{the significance level of the rule used with
#'   \code{method = "auto"}.}
#' \item{interpretation}{a short interpretation of the first imputation.}
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
#' Scholz, F. and Zhu, A. (2025) \emph{kSamples: K-Sample Rank Tests and their
#' Combinations}. R package version 1.2-12.
#' \url{https://CRAN.R-project.org/package=kSamples}
#'
#' Srivastava, M.S. and Dolatabadi, M. (2009) Multiple imputation and other
#' resampling schemes for imputing missing observations. \emph{Journal of
#' Multivariate Analysis}, \bold{100}, 1919--1937.
#'
#' @seealso \code{\link{mcar_little_test}}, \code{\link{mcar_logistic_test}},
#'   \code{\link{geewa}}, \code{\link{geewa_binary}}.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(240), ncol = 3)
#' x[41:60, 3] <- NA
#' x[61:80, 2:3] <- NA
#' mcar_homoscedasticity_test(x, method = "nonparametric")
#'
#' ## Several imputations
#' out <- mcar_homoscedasticity_test(x, method = "nonparametric",
#'                                   n_imputations = 5)
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
    n_replicates = 10000L,
    n_min = 30L,
    sig_level = 0.05,
    seed = 110L,
    maxiter = 1000L,
    tolerance = 1e-8) {
  data_name <- deparse1(substitute(object))
  imputation_supplied <- !missing(imputation)
  method <- match_mcar_method(method[1L])
  imputation <- match_mcar_imputation(imputation[1L])

  n_imputations <- check_integer_at_least(n_imputations, "n_imputations")
  min_pattern_size <- check_integer_at_least(
    min_pattern_size,
    "min_pattern_size",
    lower = 2L
  )
  n_replicates <- check_integer_at_least(n_replicates, "n_replicates")
  n_min <- check_integer_at_least(n_min, "n_min", lower = 2L)
  check_probability_open(sig_level, "sig_level")
  if (!is.null(seed) &&
      (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
       abs(seed) > .Machine$integer.max || seed != round(seed))) {
    stop(
      "'seed' must be NULL or a whole number within the integer range",
      call. = FALSE
    )
  }
  maxiter <- check_integer_at_least(maxiter, "maxiter")
  if (!is_positive_scalar(tolerance)) {
    stop("'tolerance' must be a positive finite number", call. = FALSE)
  }

  user_imputed <- !is.null(imputed_data)
  if (user_imputed && imputation_supplied) {
    stop(
      "'imputation' must not be supplied when 'imputed_data' is supplied",
      call. = FALSE
    )
  }
  if (user_imputed && n_imputations != 1L) {
    stop(
      "'n_imputations' must be 1 when 'imputed_data' is supplied",
      call. = FALSE
    )
  }

  if (inherits(object, "geer")) {
    object <- check_geer_object(object)
    x <- extract_mcar_matrix(object, data = data)
    data_name <- "repeated response from fitted geer object"
  } else {
    if (!is.null(data)) {
      stop(
        "'data' is only used when 'object' is a fitted 'geer' object",
        call. = FALSE
      )
    }
    x <- object
  }
  x <- check_mcar_jj_matrix(x)
  if (!anyNA(x)) {
    stop(
      "the Jamshidian-Jalal diagnostic requires missing values",
      call. = FALSE
    )
  }
  if (user_imputed) {
    imputed_data <- check_mcar_jj_imputed(imputed_data, x)
  }

  pattern_info <- compute_mcar_patterns(
    x,
    min_pattern_size = min_pattern_size
  )
  x_used <- check_mcar_jj_matrix(pattern_info$x)
  attr(x_used, "row.index") <- NULL
  attr(x_used, "original.nrow") <- NULL
  if (!anyNA(x_used)) {
    stop(
      paste0(
        "the Jamshidian-Jalal diagnostic requires missing values in at least ",
        "one retained pattern"
      ),
      call. = FALSE
    )
  }

  if (!user_imputed && identical(imputation, "normal") &&
      !identical(method, "hawkins")) {
    warning(
      paste0(
        "normal-theory imputation can inflate the size of the nonparametric ",
        "test for non-normal data; consider imputation = \"distribution-free\""
      ),
      call. = FALSE
    )
  }

  run_hawkins <- !identical(method, "nonparametric")
  run_nonparametric <- method %in% c("auto", "nonparametric")

  result <- run_mcar_with_seed(seed, {
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
      setup <- build_mcar_imputation_setup(
        x_used,
        imputation = imputation,
        maxiter = maxiter,
        tolerance = tolerance
      )
      completed <- lapply(
        seq_len(n_imputations),
        function(i) compute_mcar_imputation_draw(x_used, setup)
      )
    }

    neyman_nulls <- if (run_hawkins) {
      compute_mcar_neyman_nulls(
        pattern_info$group_counts,
        n_replicates = n_replicates,
        n_min = n_min
      )
    } else {
      NULL
    }

    analyses <- lapply(
      completed,
      function(completed_data) {
        hawkins <- compute_mcar_hawkins_test(
          completed = completed_data,
          group = pattern_info$group,
          group_counts = pattern_info$group_counts,
          neyman_nulls = neyman_nulls,
          test_uniformity = run_hawkins
        )
        nonparametric <- NULL
        if (run_nonparametric) {
          nonparametric <- compute_mcar_np_test(
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
    method_label <- paste0(
      "Jamshidian-Jalal ",
      "modified Hawkins MCAR screening diagnostic"
    )
    selected_test <- "hawkins"
  } else if (identical(method, "nonparametric")) {
    statistic <- c("Anderson-Darling" = nonparametric$statistic)
    parameter <- NULL
    p_value <- nonparametric$p.value
    method_label <- paste0(
      "Jamshidian-Jalal ",
      "nonparametric MCAR screening diagnostic"
    )
    selected_test <- "nonparametric"
  } else if (hawkins$p.value > sig_level) {
    statistic <- c("Fisher chi-squared" = hawkins$statistic)
    parameter <- c(df = hawkins$parameter)
    p_value <- hawkins$p.value
    method_label <- paste0(
      "Jamshidian-Jalal ",
      "automatic MCAR screening diagnostic (modified Hawkins)"
    )
    selected_test <- "hawkins"
  } else {
    statistic <- c("Anderson-Darling" = nonparametric$statistic)
    parameter <- NULL
    p_value <- nonparametric$p.value
    method_label <- paste0(
      "Jamshidian-Jalal ",
      "automatic MCAR screening diagnostic (nonparametric)"
    )
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
        if (is.null(analysis[[test]])) {
          NA_real_
        } else {
          analysis[[test]][[component]]
        }
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
        lapply(
          result$analyses,
          function(analysis) analysis$hawkins$group.p.values
        )
      )
    } else {
      NULL
    },
    nonparametric.group.statistics = if (run_nonparametric) {
      do.call(
        rbind,
        lapply(
          result$analyses,
          function(analysis) analysis$nonparametric$group.statistics
        )
      )
    } else {
      NULL
    }
  )
  if (!is.null(imputations$hawkins.group.p.values)) {
    dimnames(imputations$hawkins.group.p.values) <- list(NULL, pattern_names)
  }
  if (!is.null(imputations$nonparametric.group.statistics)) {
    dimnames(imputations$nonparametric.group.statistics) <- list(
      NULL,
      pattern_names
    )
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
    sig.level = sig_level,
    interpretation = build_mcar_interpretation(
      method = method,
      hawkins = hawkins,
      nonparametric = nonparametric,
      sig_level = sig_level
    )
  )
  class(out) <- c("mcar_homoscedasticity_test", "geer_htest", "htest")
  out
}
