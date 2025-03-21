% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wald_test.R
\name{wald_test}
\alias{wald_test}
\title{Wald Test of Nested GEE Models}
\usage{
wald_test(object0, object1, cov_type = "robust")
}
\arguments{
\item{object0}{A fitted model of the class \code{geer}.}

\item{object1}{A fitted model of the class \code{geer}.}

\item{cov_type}{character indicating whether the sandwich (robust)
covariance
matrix (\code{cov_type = "robust"}), the model-based (naive) covariance
matrix (\code{cov_type = "naive"}), the bias-corrected covariance
matrix (\code{cov_type = "bias-corrected"}) or the degrees of freedom adjusted
covariance matrix (\code{cov_type = "df-adjusted"}) should be used. By
default, the robust covariance matrix is used.}
}
\description{
Comparing two nested GEE models by carrying out a Wald test.
}
\details{
The two GEE models implied by \code{object0} and \code{object1} must be
nested.
}
\examples{
data("cerebrovascular")
fitted_model <- geewa_binary(formula = ecg ~ period * treatment,
                             id = id,
                             data = cerebrovascular,
                             link = "logit",
                             orstr = "exchangeable",
                             method = "gee")
reduced_model <- update(fitted_model, formula = ecg ~ period)
wald_test(fitted_model, reduced_model)

}
\author{
Anestis Touloumis
}
