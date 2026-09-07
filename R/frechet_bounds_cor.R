#' @title
#' Frechet Bounds for a Working Correlation Matrix
#'
#' @description
#' For a fitted \code{geer} model from \code{\link{geewa}} with a
#' \code{binomial} or \code{quasibinomial} family, a binary response and a
#' non-independence association structure, checks whether each off-diagonal
#' entry of the working correlation matrix lies within the Frechet bounds
#' implied by the fitted marginal probabilities. Results are summarized at the
#' time-pair level.
#'
#' @param object an object of class \code{geer} fitted via \code{\link{geewa}}
#'   with \code{family = binomial()} or \code{family = quasibinomial()}, a
#'   binary (0/1) response and a non-independence association structure.
#'
#' @details
#' For a pair of observations at times \eqn{j} and \eqn{k} within cluster
#' \eqn{i}, with fitted marginal probabilities \eqn{\pi_{ij}} and
#' \eqn{\pi_{ik}}, the Frechet bounds on their correlation are
#' \deqn{
#'   \ell_{ijk} = \max\!\left(
#'     -\sqrt{\frac{\pi_{ij}\pi_{ik}}{(1-\pi_{ij})(1-\pi_{ik})}},\;
#'     -\sqrt{\frac{(1-\pi_{ij})(1-\pi_{ik})}{\pi_{ij}\pi_{ik}}}
#'   \right)
#' }
#' \deqn{
#'   u_{ijk} = \min\!\left(
#'     \sqrt{\frac{\pi_{ij}(1-\pi_{ik})}{\pi_{ik}(1-\pi_{ij})}},\;
#'     \sqrt{\frac{\pi_{ik}(1-\pi_{ij})}{\pi_{ij}(1-\pi_{ik})}}
#'   \right)
#' }
#' The working correlation value \code{alpha_value} for a time pair
#' \eqn{(j, k)} is the same for all clusters and is read from the fitted
#' working correlation matrix. The bounds \eqn{\ell_{ijk}} and \eqn{u_{ijk}}
#' vary across clusters because they depend on the cluster-specific fitted
#' probabilities. The tightest bounds across clusters, \code{lower_max}
#' (maximum lower bound) and \code{upper_min} (minimum upper bound), are
#' reported in the returned data frame. The column \code{n_violated} counts
#' the number of clusters for which \code{alpha_value} falls outside
#' \eqn{[\ell_{ijk},\, u_{ijk}]}, out of the \code{n_clusters} clusters that
#' contribute an observation at both times.
#'
#' Both bounds are ratios that are undefined when a fitted probability equals
#' \eqn{0} or \eqn{1}, because the corresponding marginal distribution is
#' degenerate and its correlation with any other variable does not exist. Such
#' fitted values are therefore rejected with an error rather than propagated
#' as \code{NaN}. They typically indicate separation in the mean model.
#'
#' The bounds also assume Bernoulli marginals, so the response must be binary.
#' Grouped binomial responses supplied as a two-column matrix of successes and
#' failures are rejected, since their fitted values are Binomial proportions
#' rather than Bernoulli means.
#'
#' Clusters of size one contribute no time pair and are skipped. A time pair
#' observed in no cluster produces no row.
#'
#' @return
#' A data frame with one row per unique time pair \eqn{(j, k)} and columns:
#' \describe{
#'   \item{\code{alpha_name}}{label of the form \code{alpha_j.k} identifying the
#'     working correlation entry for this time pair.}
#'   \item{\code{alpha_value}}{working correlation value for the time pair.}
#'   \item{\code{lower_max}}{maximum Frechet lower bound across clusters,
#'     giving the tightest lower admissibility constraint.}
#'   \item{\code{upper_min}}{minimum Frechet upper bound across clusters,
#'     giving the tightest upper admissibility constraint.}
#'   \item{\code{n_clusters}}{number of clusters contributing an observation at
#'     both times, and hence the number of bound pairs summarized in
#'     \code{lower_max} and \code{upper_min}.}
#'   \item{\code{n_violated}}{number of clusters for which \code{alpha_value}
#'     falls outside the cluster-specific Frechet bounds.}
#' }
#'
#' @seealso
#' \code{\link{geewa}}.
#'
#' @examples
#' data("cholecystectomy", package = "geer")
#'
#' fit <- geewa(
#'   formula = pain ~ treatment + gender + age + I(time > 4),
#'   family = binomial(link = "logit"),
#'   data = cholecystectomy,
#'   id = id,
#'   repeated = time,
#'   corstr = "unstructured",
#'   method = "gee"
#' )
#' frechet_bounds_cor(fit)
#'
#' @export
frechet_bounds_cor <- function(object) {
  object <- check_geer_object(object)
  if (!identical(object$fit_function, "geewa")) {
    stop(
      "'object' must be fitted by 'geewa': 'geewa_binary' parameterizes the ",
      "within-cluster association through marginalized odds ratios rather than ",
      "correlations, so the Frechet bounds on a working correlation matrix do ",
      "not apply",
      call. = FALSE
    )
  }
  if (!(object$family$family %in% c("binomial", "quasibinomial"))) {
    stop(
      "'object' must be a binomial fit: got '", object$family$family, "'",
      call. = FALSE
    )
  }
  if (identical(object$association_structure, "independence")) {
    stop(
      "'object' must not have an 'independence' association structure",
      call. = FALSE
    )
  }
  y <- object$y
  if (anyNA(y) || !all(y == 0 | y == 1)) {
    stop(
      "'object' must have a binary (0/1) response: the Frechet bounds assume ",
      "Bernoulli marginals, so a grouped binomial response supplied as a ",
      "two-column matrix of successes and failures is not supported",
      call. = FALSE
    )
  }
  mu <- object$fitted.values
  degenerate <- !is.finite(mu) | mu <= 0 | mu >= 1
  if (any(degenerate)) {
    stop(
      sprintf(
        paste0(
          "the Frechet bounds are undefined: %d of %d fitted probabilities ",
          "are not strictly between 0 and 1, so the corresponding marginal ",
          "distributions are degenerate. This usually indicates separation in ",
          "the mean model"
        ),
        sum(degenerate),
        length(mu)
      ),
      call. = FALSE
    )
  }
  id <- object$id
  repeated <- as.integer(object$repeated)
  time_max <- max(repeated)
  cor_mat <- get_correlation_matrix(
    object$association_structure,
    object$alpha,
    time_max
  )
  cluster_index <- split(seq_along(id), id)
  cluster_index <- cluster_index[lengths(cluster_index) >= 2L]
  if (length(cluster_index) == 0L) {
    return(
      data.frame(
        alpha_name = character(0),
        alpha_value = numeric(0),
        lower_max = numeric(0),
        upper_min = numeric(0),
        n_clusters = integer(0),
        n_violated = integer(0),
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    )
  }
  pair_index <- lapply(
    cluster_index,
    function(idx) {
      combinations <- utils::combn(length(idx), 2L)
      rbind(idx[combinations[1L, ]], idx[combinations[2L, ]])
    }
  )
  pair_index <- do.call(cbind, pair_index)
  first <- pair_index[1L, ]
  second <- pair_index[2L, ]
  p <- mu[first]
  q <- mu[second]
  lower <- pmax(
    -sqrt(p * q / ((1 - p) * (1 - q))),
    -sqrt((1 - p) * (1 - q) / (p * q))
  )
  upper <- pmin(
    sqrt(p * (1 - q) / (q * (1 - p))),
    sqrt(q * (1 - p) / (p * (1 - q)))
  )
  time_j <- repeated[first]
  time_k <- repeated[second]
  alpha_value <- cor_mat[cbind(time_j, time_k)]
  violated <- alpha_value < lower | alpha_value > upper
  pair_key <- factor((time_j - 1L) * time_max + time_k)
  key_values <- as.integer(levels(pair_key))
  out_j <- (key_values - 1L) %/% time_max + 1L
  out_k <- (key_values - 1L) %% time_max + 1L
  data.frame(
    alpha_name = paste0("alpha_", out_j, ".", out_k),
    alpha_value = cor_mat[cbind(out_j, out_k)],
    lower_max = as.numeric(tapply(lower, pair_key, max)),
    upper_min = as.numeric(tapply(upper, pair_key, min)),
    n_clusters = as.integer(tapply(violated, pair_key, length)),
    n_violated = as.integer(tapply(violated, pair_key, sum)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
