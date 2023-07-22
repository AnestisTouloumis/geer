#ifndef COVARIANCEMATRICES_H
#define COVARIANCEMATRICES_H

Rcpp::List get_covariance_matrices(const arma::vec & y_vector,
                                   const arma::mat & model_matrix,
                                   const arma::vec & id_vector,
                                   const arma::vec & repeated_vector,
                                   const char * link,
                                   const char * family,
                                   const arma::vec & mu_vector,
                                   const arma::vec & eta_vector,
                                   const char * correlation_structure,
                                   const arma::vec & alpha_vector,
                                   const double & phi);

#endif
