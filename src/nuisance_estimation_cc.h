#ifndef GEER_NUISANCE_ESTIMATION_CC_H
#define GEER_NUISANCE_ESTIMATION_CC_H


#include <RcppArmadillo.h>
#include <vector>
#include "cluster_utils.h"
#include "family_codes.h"

arma::vec get_pearson_residuals(const char* family,
                                const arma::vec& y_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& weights_vector);

arma::vec get_pearson_residuals(FamilyCode family_code,
                                const arma::vec& y_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& weights_vector);
double get_phi_hat_cc(const arma::vec& pearson_residuals_vector,
                   const int params_no);
arma::vec get_alpha_hat_cc(const char* correlation_structure,
                        const arma::vec& pearson_residuals_vector,
                        const std::vector<Cluster>& clusters,
                        const arma::vec& repeated_vector,
                        const double phi,
                        const int params_no,
                        const int mdependence);

#endif
