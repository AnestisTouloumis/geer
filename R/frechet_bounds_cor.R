#' @title
#' \enc{Fréchet}{Frechet} Bounds for a Working Correlation Matrix
#'
#' @description
#' Checks whether each working correlation of a fitted binary GEE model lies
#' within the \enc{Fréchet}{Frechet} bounds implied by the fitted marginal
#' probabilities.
#'
#' @param object a fitted model object of class \code{"geer"}, obtained from
#'   \code{\link{geewa}} with a \code{binomial} or \code{quasibinomial} family,
#'   a binary (0/1) response, and a working correlation structure other than
#'   \code{"independence"}.
#'
#' @details
#' For observations at times \eqn{j} and \eqn{k} of cluster \eqn{i}, with
#' fitted probabilities \eqn{\pi_{ij}}{pi_ij} and \eqn{\pi_{ik}}{pi_ik}, the
#' \enc{Fréchet}{Frechet} bounds on their correlation are
#' \deqn{\ell_{ijk} = \max\left(
#'   -\sqrt{\frac{\pi_{ij}\pi_{ik}}{(1-\pi_{ij})(1-\pi_{ik})}},
#'   -\sqrt{\frac{(1-\pi_{ij})(1-\pi_{ik})}{\pi_{ij}\pi_{ik}}}\right),}{l_ijk = max(-sqrt(pi_ij pi_ik/((1 - pi_ij)(1 - pi_ik))), -sqrt((1 - pi_ij)(1 - pi_ik)/(pi_ij pi_ik))),}
#' \deqn{u_{ijk} = \min\left(
#'   \sqrt{\frac{\pi_{ij}(1-\pi_{ik})}{\pi_{ik}(1-\pi_{ij})}},
#'   \sqrt{\frac{\pi_{ik}(1-\pi_{ij})}{\pi_{ij}(1-\pi_{ik})}}\right).}{u_ijk = min(sqrt(pi_ij (1 - pi_ik)/(pi_ik (1 - pi_ij))), sqrt(pi_ik (1 - pi_ij)/(pi_ij (1 - pi_ik)))).}
#' The working correlation of a time pair is common to all clusters, whereas
#' the bounds vary with the fitted probabilities. The bounds are undefined for
#' fitted probabilities of 0 or 1, which usually indicate separation, and such
#' fits give an error. Grouped binomial responses are not allowed. Clusters of
#' size one are skipped, and time pairs observed in no cluster are omitted.
#'
#' @return
#' A data frame with one row per time pair and columns:
#' \item{alpha_name}{the label \code{alpha_j.k} of the working correlation.}
#' \item{alpha_value}{the working correlation.}
#' \item{lower_max}{the largest lower bound over the clusters.}
#' \item{upper_min}{the smallest upper bound over the clusters.}
#' \item{n_clusters}{the number of clusters observed at both times.}
#' \item{n_violated}{the number of clusters whose bounds do not contain
#'   \code{alpha_value}.}
#'
#' @seealso \code{\link{geewa}}.
#'
#' @examples
#' data("cholecystectomy", package = "geer")
#' fit <- geewa(
#'   formula = pain ~ treatment + gender + age + I(time > 4),
#'   family = binomial(link = "logit"),
#'   data = cholecystectomy,
#'   id = id,
#'   repeated = time,
#'   corstr = "unstructured"
#' )
#' frechet_bounds_cor(fit)
#'
#' @export
frechet_bounds_cor <- function(object) {
  object <- check_geer_object(object)
  if (!identical(get_geer_fit_function(object), "geewa")) {
    stop(
      paste0(
        "'object' must be a geewa() fit: the Frechet bounds apply to working ",
        "correlations, which geewa_binary() fits do not use"
      ),
      call. = FALSE
    )
  }
  if (!(object$family$family %in% c("binomial", "quasibinomial"))) {
    stop(
      sprintf(
        "'object' must be a binomial fit, not '%s'",
        object$family$family
      ),
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
  if (!is_binary_response(y)) {
    stop(
      paste0(
        "'object' must have a binary (0/1) response; grouped binomial ",
        "responses are not supported"
      ),
      call. = FALSE
    )
  }
  mu <- object$fitted.values
  degenerate <- is_degenerate_probability(mu)
  if (any(degenerate)) {
    stop(
      sprintf(
        paste0(
          "the Frechet bounds are undefined: %d of %d fitted probabilities ",
          "are 0 or 1, which suggests separation"
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
  cor_mat <- get_correlation_matrix_cc(
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
