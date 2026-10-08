#ifndef GEER_WORKING_COVARIANCE_OR_H
#define GEER_WORKING_COVARIANCE_OR_H

// Working covariance matrix V_i of the odds-ratio GEE for binary data, its
// derivatives with respect to the means, and the matrices G_i that enter them.

#include <RcppArmadillo.h>

arma::vec get_subject_specific_odds_ratios(const arma::vec& repeated_vector_i,
                                           const arma::uword cluster_size_max,
                                           const arma::vec& odds_ratios_vector);
arma::mat get_v_matrix_or(const arma::vec& mu_vector,
                          const arma::vec& odds_ratios_vector,
                          const arma::vec& weights_vector);
arma::mat get_g_matrix(const arma::vec& mu_vector,
                       const arma::vec& odds_ratios_vector);
arma::mat get_g_matrix_mu(const arma::vec& mu_vector,
                          const arma::vec& odds_ratios_vector);
arma::mat get_v_matrix_mu_or(const arma::vec& mu_vector,
                             const arma::vec& odds_ratios_vector,
                             const arma::vec& weights_vector);

#endif
