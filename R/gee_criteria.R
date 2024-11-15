#' Variable and Covariance Selection Criteria
#'
#' Reports commonly used criteria for variable selection
#' and for selecting the "working" association structure for one or several
#' fitted models from the \code{geer} package.
#'
#'
#' The Quasi Information Criterion (QIC), the Correlation Information Criterion
#' (CIC) and the Rotnitzky and Jewell Criterion (RJC) are used for selecting the
#' best association structure. The QICu criterion is used for selecting the best
#' subset of covariates. When choosing among GEE models with different association
#' structures but with the same subset of covariates, the model with the smallest
#' value of QIC, CIC or RJC should be preffered. When choosing between GEE models
#' with different number of covariates, the model with the smallest QICu value
#' should be preferred.
#'
#' @return A vector or matrix with the QIC, QICu, CIC, RJC and the number of
#' regression parameters (including intercepts).
#'
#' @param object an object of the class \code{geer}.
#' @param ... optionally more objects of the class \code{geer}.
#' @param cov_type character indicating whether the sandwich (robust)
#' covariance
#' matrix (\code{cov_type = "robust"}), the bias-corrected covariance
#' matrix (\code{cov_type = "bias-corrected"}) or the degrees of freedom adjusted
#' covariance matrix (\code{cov_type = "df-adjusted"}) should be used in computing
#' the estimated covariance matrix. By default, the robust covariance matrix is
#' used.
#'
#' @author Anestis Touloumis
#' @references Hin, L.Y. and Wang, Y.G. (2009) Working correlation structure
#' identification in generalized estimating equations. \emph{Statistics in
#' Medicine} \bold{28}, 642--658.
#'
#' Pan, W. (2001) Akaike's information criterion in generalized
#' estimating equations. \emph{Biometrics} \bold{57}, 120--125.
#'
#' Rotnitzky, A. and Jewell, N.P. (1990) Hypothesis testing of regression
#' parameters in semiparametric generalized linear models for cluster correlated
#' data. \emph{Biometrika} \bold{77}, 485--497.
#'
#' @examples
#' data("cerebrovascular")
#' fitted_model <- geewa_binary(formula = ecg ~ period * treatment,
#'                              id = id,
#'                              data = cerebrovascular,
#'                              link = "logit",
#'                              or_structure = "exchangeable",
#'                              method = "gee")
#' gee_criteria(fitted_model)
#'
#' @export
gee_criteria <- function(object, ..., cov_type = "robust") {
  UseMethod("gee_criteria")
}

#' @aliases gee_criteria gee_criteria.geer
#'
#' @method gee_criteria geer
#' @export

gee_criteria.geer <- function(object, ..., cov_type = "robust") {
  if (length(list(...))) {
    results <- lapply(list(object, ...), compute_criteria, cov_type)
    #check <- sapply(list(object, ...), function(x) {
     # length(x$y)
    #})
   # if (any(check != check[1]))
     # warning("models are not all fitted to the same number of observations")
    res <- do.call("rbind", results)
    Call <- match.call()
    if("cov_type" %in% names(col))
    #row.names(res) <- as.character(Call[-(length(list(...)))])
    Call
  } else {
    compute_criteria(object, cov_type)
  }
}


compute_criteria <- function(object, cov_type = "robust") {
  mu <- object$fitted.values
  y  <- object$y
  marginal_distribution <- object$family$family
  quasi_likelihood <-
    switch(marginal_distribution,
           gaussian = sum(((y - mu)^2)/-2),
           binomial = sum(y * log(mu) + (1 - y) * log(1 - mu)),
           poisson  = sum((y * log(mu)) - mu),
           Gamma    = sum(-y/(mu - log(mu))),
           stop("Error: distribution not recognized")
    )/object$phi
  s_vector <- y - mu
  if (any(names(object$call) == "correlation_structure")) {
    independence_model <- update(object,
                                 correlation_structure = "independence")
  } else {
    independence_model <- update(object,
                                 or_structure = "independence")
  }
  independence_naive_covariance <- vcov(independence_model,
                                        cov_type = "naive")
  naive_covariance <- vcov(object,
                           cov_type = "naive")
  robust_covariance <- vcov(object,
                            cov_type = cov_type)
  p <- length(object$coeff)
  qic_u <- 2 * (p - quasi_likelihood)
  cic <- sum(solve(independence_naive_covariance) * robust_covariance)
  qic <- 2 * (cic - quasi_likelihood)
  q_matrix <- robust_covariance %*% solve(naive_covariance)
  c1 <- sum(diag(q_matrix)) / p
  c2 <- sum(q_matrix * t(q_matrix)) / p
  rjc <- sqrt(((1 - c1) ^ 2) + ((1 - c2) ^ 2))
  ans <- c(qic, cic, rjc, qic_u, p)
  names(ans) <- c("QIC", "CIC", "RJC", "QICu", "Parameters")
  ans
}
