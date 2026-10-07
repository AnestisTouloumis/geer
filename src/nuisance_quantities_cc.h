#ifndef NUISANCE_QUANTITIES_CC_H
#define NUISANCE_QUANTITIES_CC_H

#include <RcppArmadillo.h>
#include "family_codes.h"

// char* overload (R-facing)
arma::vec get_pearson_residuals(const char* family,
                                const arma::vec& y_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& weights_vector);

// FamilyCode overload (hot-path)
arma::vec get_pearson_residuals(FamilyCode fc,
                                const arma::vec& y_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& weights_vector);
double get_phi_hat(const arma::vec& pearson_residuals_vector,
                   const int& params_no);
arma::vec get_alpha_hat(const char* correlation_structure,
                        const arma::vec& pearson_residuals_vector,
                        const arma::vec& id_vector,
                        const arma::vec& repeated_vector,
                        const double& phi,
                        const int& params_no,
                        const int& mdependence);
arma::mat get_correlation_matrix(const char* correlation_structure,
                                 const arma::vec& alpha_vector,
                                 const arma::uword dimension);

arma::mat get_v_matrix_cc(FamilyCode fc,
                          const arma::vec& mu_vector,
                          const arma::vec& repeated_vector,
                          const double& phi,
                          const arma::mat& cor_matrix,
                          const arma::vec& weights_vector);

arma::mat get_v_matrix_cc(const char* family,
                          const arma::vec& mu_vector,
                          const arma::vec& repeated_vector,
                          const double& phi,
                          const arma::mat& cor_matrix,
                          const arma::vec& weights_vector);


#endif
