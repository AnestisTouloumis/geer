#ifndef VARIANCE_FUNCTIONS_H
#define VARIANCE_FUNCTIONS_H

#include <RcppArmadillo.h>
#include "family_utils.h"

arma::vec variance(FamilyCode fc,    const arma::vec& mu_vector);
arma::vec variancemu(FamilyCode fc,  const arma::vec& mu_vector);
arma::vec variancemu2(FamilyCode fc, const arma::vec& mu_vector);

arma::vec variance(const char* family,    const arma::vec& mu_vector);
arma::vec variancemu(const char* family,  const arma::vec& mu_vector);
arma::vec variancemu2(const char* family, const arma::vec& mu_vector);

#endif
