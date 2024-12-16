#ifndef NUISANCEQUANTITIES_H
#define NUISANCEQUANTITIES_H

double get_phi_hat(const arma::vec & pearson_residuals_vector,
                   const int & params_no);

double alpha_hat_exchangeable(const arma::vec & pearson_residuals_vector,
                              const arma::vec & id_vector,
                              const double & phi,
                              const int & params_no);
double alpha_hat_ar1(const arma::vec & pearson_residuals_vector,
                     const arma::vec & id_vector,
                     const arma::vec & repeated_vector,
                     const double & phi,
                     const int & params_no);
arma::vec alpha_hat_unstructured(const arma::vec & pearson_residuals_vector,
                                 const arma::vec & id_vector,
                                 const arma::vec & repeated_vector,
                                 const double & phi,
                                 const int & params_no);
arma::vec alpha_hat_mdependent(const arma::vec & pearson_residuals_vector,
                               const arma::vec & id_vector,
                               const arma::vec & repeated_vector,
                               const double & phi,
                               const int & params_no,
                               const int & mdependence);
arma::vec alpha_hat_toeplitz(const arma::vec & pearson_residuals_vector,
                             const arma::vec & id_vector,
                             const arma::vec & repeated_vector,
                             const double & phi,
                             const int & params_no);
arma::vec get_alpha_hat(const char * correlation_structure,
                        const arma::vec & pearson_residuals_vector,
                        const arma::vec & id_vector,
                        const arma::vec & repeated_vector,
                        const double & phi,
                        const int & params_no,
                        const int & mdependence);

arma::mat correlation_independence(const int & dimension);
arma::mat correlation_exchangeable(const arma::vec & alpha_vector,
                                   const int & dimension);
arma::mat correlation_exchangeable_inverse(const arma::vec & alpha_vector,
                                           const int & dimension);
arma::mat correlation_ar1(const arma::vec & alpha_vector,
                          const int & dimension);
arma::mat correlation_ar1_inverse(const arma::vec & alpha_vector,
                                  const int & dimension);
arma::mat correlation_mdependent(const arma::vec & alpha_vector,
                                 const int & dimension);
arma::mat correlation_toeplitz(const arma::vec & alpha_vector);
arma::mat correlation_unstructured(const arma::vec & alpha_vector,
                                   const int & dimension);

arma::mat get_correlation_matrix(const char * correlation_structure,
                                 const arma::vec & alpha_vector,
                                 const int & dimension);
arma::mat get_correlation_matrix_inverse(const char * correlation_structure,
                                         const arma::vec & alpha_vector,
                                         const int & dimension);

arma::vec get_pearson_residuals(const char * family,
                                const arma::vec & y_vector,
                                const arma::vec & mu_vector);

arma::mat get_weight_matrix(const char * family,
                            const arma::vec & mu_vector,
                            const arma::vec & repeated_vector,
                            const double & phi,
                            const arma::mat & cor_matrix_inverse,
                            const arma::vec & weights_vector);

#endif
