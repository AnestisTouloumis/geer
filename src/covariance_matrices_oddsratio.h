#ifndef COVARIANCEMATRICESODDSRATIO_H
#define COVARIANCEMATRICESODDSRATIO_H

Rcpp::List get_covariance_matrices_or(const arma::vec & y_vector,
                                      const arma::mat & model_matrix,
                                      const arma::vec & id_vector,
                                      const arma::vec & repeated_vector,
                                      const char * link,
                                      const arma::vec & mu_vector,
                                      const arma::vec & eta_vector,
                                      const arma::vec & alpha_vector);

#endif
