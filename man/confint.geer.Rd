% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geer_class.R
\name{confint.geer}
\alias{confint.geer}
\alias{confint}
\title{Confidence Intervals for Model Parameters from a \code{geer} Object}
\usage{
\method{confint}{geer}(object, parm, level = 0.95, cov_type = "robust", ...)
}
\arguments{
\item{object}{an object representing a model of the class \code{geer}.}

\item{parm}{a specification of which parameters are to be given
    confidence intervals, either a vector of numbers or a vector of
    names.  If missing, all parameters are considered.}

\item{level}{the confidence level required.}

\item{cov_type}{character indicating the covariance matrix estimator
of the regression parameters. Options include the sandwich or
robust estimator (\code{"robust"}), the bias-corrected estimator
(\code{"bias-corrected"}), the degrees of freedom adjusted estimator
(\code{"df-adjusted"}) and the model-based or naive estimator
(\code{"naive"}). By default, the robust covariance matrix
estimator is used.}

\item{...}{additional argument(s) for methods.}
}
\value{
A matrix (or vector) with columns giving lower and upper confidence
  limits for each parameter. These will be labelled as (1-level)/2 and
  1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
}
\description{
Confidence intervals for one or more parameters in a \code{geer} object.
}
\details{
The references in \code{\link{vcov}} include the formulae for the covariance
matrix implied by \code{cov_type}.
}
\examples{
data("leprosy")
fit <- geewa(formula = bacilli ~ factor(period) + factor(period):treatment,
             family = poisson(link = "log"), id = id, data = leprosy)
confint(fit)

data("cerebrovascular")
fit <- geewa_binary(formula = ecg ~ treatment + factor(period), link = "logit",
                    id = id, data = cerebrovascular)
confint(fit)

}
