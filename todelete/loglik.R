object <- object2
object <- update(object2, family = gaussian())

fam <- object$family$family
p <- length(object$coef)
dev <- residuals(object, "deviance")
## allow for estimated dispersion
if(fam %in% c("gaussian", "Gamma", "inverse.gaussian")) p <- p + 1
loglik <- -object$family$aic(object$y,
                         object$obs_no,
                         fitted(object),
                         object$weights,
                         dev)/2
loglik
compute_quasi_loglikelihood(object)
