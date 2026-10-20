#ifndef GEER_INFLUENCE_MEASURES_H
#define GEER_INFLUENCE_MEASURES_H

#include <RcppArmadillo.h>


Rcpp::List get_influence_quantities_cc(const arma::vec& y_vector,
                                       const arma::mat& model_matrix,
                                       const arma::vec& id_vector,
                                       const arma::vec& repeated_vector,
                                       const arma::vec& weights_vector,
                                       const char* link,
                                       const char* family,
                                       const arma::vec& mu_vector,
                                       const arma::vec& eta_vector,
                                       const char* correlation_structure,
                                       const arma::vec& alpha_vector,
                                       const double phi);
Rcpp::List get_influence_quantities_or(const arma::vec& y_vector,
                                       const arma::mat& model_matrix,
                                       const arma::vec& id_vector,
                                       const arma::vec& repeated_vector,
                                       const arma::vec& weights_vector,
                                       const char* link,
                                       const arma::vec& mu_vector,
                                       const arma::vec& eta_vector,
                                       const arma::vec& alpha_vector);


#endif
