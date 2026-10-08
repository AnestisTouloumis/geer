#ifndef GEER_WORKING_COVARIANCE_CC_H
#define GEER_WORKING_COVARIANCE_CC_H

// Working correlation matrix R(alpha) and working covariance matrix V_i of the
// correlation-structure GEE.

#include <RcppArmadillo.h>
#include "family_codes.h"

arma::mat get_correlation_matrix(const char* correlation_structure,
                                 const arma::vec& alpha_vector,
                                 const arma::uword dimension);

arma::mat get_v_matrix_cc(FamilyCode family_code,
                          const arma::vec& mu_vector,
                          const arma::vec& repeated_vector,
                          const double phi,
                          const arma::mat& cor_matrix,
                          const arma::vec& weights_vector);

arma::mat get_v_matrix_cc(const char* family,
                          const arma::vec& mu_vector,
                          const arma::vec& repeated_vector,
                          const double phi,
                          const arma::mat& cor_matrix,
                          const arma::vec& weights_vector);

#endif
