% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geer_class.R
\name{model.matrix.geer}
\alias{model.matrix.geer}
\alias{model.matrix}
\title{Construct Design Matrices from a \code{geer} Object}
\usage{
\method{model.matrix}{geer}(object, ...)
}
\arguments{
\item{object}{an object representing a model of the class \code{geer}.}

\item{...}{additional argument(s) for methods.}
}
\value{
The design matrix for the marginal regression model with the specified formula and
data.

There is an attribute \code{"assign"}, an integer vector with an entry for
each column in the matrix giving the term in the formula which gave rise to
the column. Value 0 corresponds to the intercept (if any), and positive
values to terms in the order given by the \code{term.labels} attribute of the
terms structure corresponding to \code{object}.

If there are any factors in terms in the model, there is an attribute
\code{"contrasts"}, a named list with an entry for each factor. This
specifies the contrasts that would be used in terms in which the factor is
coded by contrasts (in some terms dummy coding may be used), either as a
character vector naming a function or as a numeric matrix.
}
\description{
Creates a design or model matrix from a \code{geer} object, e.g., by expanding
factors to a set of dummy variables (depending on the contrasts) and expanding
interactions similarly.
}
\details{
\code{model.matrix} creates a design matrix from the description given in
\code{terms(object)}, using the data in \code{object$data}.

In an interaction term, the variable whose levels vary fastest is the first
one to appear in the formula (and not in the term), so in
\code{~ a + b + b:a} the interaction will have a varying fastest.

By convention, if the response variable also appears on the right-hand side
of the formula it is dropped (with a warning), although interactions
involving the term are retained.
}
\examples{
data("leprosy")
fit <- geewa(formula = bacilli ~ factor(period) + factor(period):treatment,
             family = poisson(link = "log"), id = id, data = leprosy)
model.matrix(fit)

data("cerebrovascular")
fit <- geewa_binary(formula = ecg ~ treatment + factor(period), link = "logit",
                    id = id, data = cerebrovascular)
model.matrix(fit)

}
