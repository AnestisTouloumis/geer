#' Fréchet Bounds for a Working Correlation Matrix
#'
#' For a fitted \code{geer} model from \code{\link{geewa}} with a
#' \code{binomial} family and a non-independence association structure,
#' checks whether each off-diagonal entry of the working correlation matrix
#' lies within the Fréchet bounds implied by the fitted marginal probabilities.
#' Results are summarised at the time-pair level.
#'
#' @param object an object of class \code{geer} fitted via \code{\link{geewa}}
#'   with \code{family = binomial()} and a non-independence association
#'   structure.
#'
#' @details
#' For a pair of observations at times \eqn{j} and \eqn{k} within cluster
#' \eqn{i}, with fitted marginal probabilities \eqn{\pi_{ij}} and
#' \eqn{\pi_{ik}}, the Fréchet bounds on their correlation are
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
#' The working correlation value \code{cor} for a time pair \eqn{(j, k)} is
#' the same for all clusters and is read from the fitted working correlation
#' matrix. The bounds \eqn{\ell_{ijk}} and \eqn{u_{ijk}} vary across clusters
#' because they depend on the cluster-specific fitted probabilities. The tightest
#' bounds across clusters, lower_max (maximum lower bound) and upper_min
#' (minimum upper bound), are reported in the returned data frame. The column
#' \code{n_violated} counts the number of clusters for which \code{cor} falls outside
#' \eqn{(\ell_{ijk},\, u_{ijk})}.
#'
#'
#' @return a data frame with one row per unique time pair \eqn{(j, k)} and
#'   columns:
#' \describe{
#'   \item{\code{alpha_name}}{label of the form \code{alpha_j.k} identifying the
#'     working correlation entry for this time pair.}
#'   \item{\code{alpha_value}}{working correlation value for the time pair.}
#'   \item{\code{lower_max}}{maximum Fréchet lower bound across clusters,
#'     giving the tightest lower admissibility constraint.}
#'   \item{\code{upper_min}}{minimum Fréchet upper bound across clusters,
#'     giving the tightest upper admissibility constraint.}
#'   \item{\code{n_violated}}{number of clusters for which \code{alpha_value}
#'     falls outside the cluster-specific Fréchet bounds.}
#' }
#'
#' @seealso \code{\link{geewa}}.
#'
#' @export
frechet_bounds_cor <- function(object) {
  if (!inherits(object, "geer")) {
    stop("'object' must be of class \"geer\".", call. = FALSE)
  }
  if (!grepl("binomial", object$family$family, ignore.case = TRUE)) {
    stop(
      "frechet_bounds_cor: family must be \"binomial\"; ",
      "got \"", object$family$family, "\".",
      call. = FALSE
    )
  }
  if (object$association_structure == "independence") {
    stop(
      "frechet_bounds_cor: association structure must not be \"independence\".",
      call. = FALSE
    )
  }
  mu       <- object$fitted.values
  id       <- object$id
  repeated <- object$repeated
  time_max <- max(repeated)
  cor_mat  <- get_correlation_matrix(
    object$association_structure,
    object$alpha,
    time_max
  )
  uid  <- unique(id)
  n_cl <- length(uid)
  pair_lower    <- list()
  pair_upper    <- list()
  pair_violated <- list()
  for (i in seq_len(n_cl)) {
    idx  <- which(id == uid[i])
    mu_i <- mu[idx]
    re_i <- repeated[idx]
    n_i  <- length(idx)
    if (n_i == 1L) next
    for (j1 in seq_len(n_i - 1L)) {
      for (j2 in (j1 + 1L):n_i) {
        p   <- mu_i[j1]
        q   <- mu_i[j2]
        tj  <- re_i[j1]
        tk  <- re_i[j2]
        key <- paste0(tj, ":", tk)
        lo  <- max(
          -sqrt(p * q / ((1 - p) * (1 - q))),
          -sqrt((1 - p) * (1 - q) / (p * q))
        )
        up  <- min(
          sqrt(p * (1 - q) / (q * (1 - p))),
          sqrt(q * (1 - p) / (p * (1 - q)))
        )
        cv  <- cor_mat[tj, tk]
        pair_lower[[key]]    <- c(pair_lower[[key]],    lo)
        pair_upper[[key]]    <- c(pair_upper[[key]],    up)
        pair_violated[[key]] <- c(pair_violated[[key]], cv < lo || cv > up)
      }
    }
  }
  keys   <- names(pair_lower)
  n_pairs <- length(keys)
  time_j     <- integer(n_pairs)
  time_k     <- integer(n_pairs)
  lower_min  <- numeric(n_pairs)
  lower_max  <- numeric(n_pairs)
  upper_min  <- numeric(n_pairs)
  upper_max  <- numeric(n_pairs)
  cor_vec    <- numeric(n_pairs)
  n_violated <- integer(n_pairs)
  for (s in seq_len(n_pairs)) {
    parts      <- as.integer(strsplit(keys[s], ":", fixed = TRUE)[[1L]])
    time_j[s]     <- parts[1L]
    time_k[s]     <- parts[2L]
    lower_min[s]  <- min(pair_lower[[keys[s]]])
    lower_max[s]  <- max(pair_lower[[keys[s]]])
    upper_min[s]  <- min(pair_upper[[keys[s]]])
    upper_max[s]  <- max(pair_upper[[keys[s]]])
    cor_vec[s]    <- cor_mat[parts[1L], parts[2L]]
    n_violated[s] <- sum(pair_violated[[keys[s]]])
  }
  corstr <- object$association_structure
  alpha_name <- paste0("alpha_", time_j, ".", time_k)
  ord <- order(time_j, time_k)
  data.frame(
    alpha_name  = alpha_name[ord],
    alpha_value = cor_vec[ord],
    lower_max   = lower_max[ord],
    upper_min   = upper_min[ord],
    n_violated  = n_violated[ord],
    row.names   = NULL
  )
}
