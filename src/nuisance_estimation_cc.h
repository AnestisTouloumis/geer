#ifndef GEER_NUISANCE_ESTIMATION_CC_H
#define GEER_NUISANCE_ESTIMATION_CC_H

// Estimation of the nuisance quantities (dispersion and correlation
// parameters) of the correlation-structure GEE: Pearson residuals, phi and
// alpha. The working covariance matrices built from them are declared in
// working_covariance_cc.h.

#include <RcppArmadillo.h>
#include <vector>
#include "cluster_utils.h"
#include "family_codes.h"

// char* overload (R-facing)
arma::vec get_pearson_residuals(const char* family,
                                const arma::vec& y_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& weights_vector);

// FamilyCode overload (hot-path)
arma::vec get_pearson_residuals(FamilyCode family_code,
                                const arma::vec& y_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& weights_vector);
double get_phi_hat(const arma::vec& pearson_residuals_vector,
                   const int params_no);
arma::vec get_alpha_hat(const char* correlation_structure,
                        const arma::vec& pearson_residuals_vector,
                        const std::vector<Cluster>& clusters,
                        const arma::vec& repeated_vector,
                        const double phi,
                        const int params_no,
                        const int mdependence);

#endif
