% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geer_class.R
\name{coef.geer}
\alias{coef.geer}
\alias{coef}
\alias{coefficients}
\title{Extract Model Coefficients from a \code{geer} Object}
\usage{
\method{coef}{geer}(object, ...)
}
\arguments{
\item{object}{an object representing a model of the class \code{geer}.}

\item{...}{additional argument(s) for methods.}
}
\value{
A named numeric vector with the coefficients extracted from \code{object}.
}
\description{
Extracts model coefficients from a \code{geer} object. \code{coefficients} is
an alias for it.
}
\examples{
data("leprosy")
fit <- geewa(formula = bacilli ~ factor(period) + factor(period):treatment,
             family = poisson(link = "log"), id = id, data = leprosy)
coef(fit)

data("cerebrovascular")
fit <- geewa_binary(formula = ecg ~ treatment + factor(period), link = "logit",
                    id = id, data = cerebrovascular)
coef(fit)

}
