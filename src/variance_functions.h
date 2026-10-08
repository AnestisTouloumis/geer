#ifndef GEER_VARIANCE_FUNCTIONS_H
#define GEER_VARIANCE_FUNCTIONS_H

#include <RcppArmadillo.h>
#include "family_codes.h"

namespace geer {

arma::vec variance_function(FamilyCode family_code, const arma::vec& mu_vector);
arma::vec variance_derivative_1(FamilyCode family_code, const arma::vec& mu_vector);
arma::vec variance_derivative_2(FamilyCode family_code, const arma::vec& mu_vector);

arma::vec variance_function(const char* family, const arma::vec& mu_vector);

}  // namespace geer

#endif
