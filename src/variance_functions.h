#ifndef VARIANCE_FUNCTIONS_H
#define VARIANCE_FUNCTIONS_H

#include <RcppArmadillo.h>

arma::vec variance(const char* family,
                   const arma::vec& mu_vector);

arma::vec variancemu(const char* family,
                     const arma::vec& mu_vector);

arma::vec variancemu2(const char* family,
                      const arma::vec& mu_vector);

#endif
