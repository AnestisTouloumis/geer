#ifndef UPDATEBETAODDSRATIO_H
#define UPDATEBETAODDSRATIO_H


arma::vec update_beta_gee_or(const arma::vec & y_vector,
                             const arma::mat & model_matrix,
                             const arma::vec & id_vector,
                             const arma::vec & repeated_vector,
                             const char * link,
                             const arma::vec & beta_vector,
                             const arma::vec & mu_vector,
                             const arma::vec & eta_vector,
                             const arma::vec & alpha_vector,
                             const arma::vec & weights_vector);

arma::vec update_beta_naive_or(const arma::vec & y_vector,
                               const arma::mat & model_matrix,
                               const arma::vec & id_vector,
                               const arma::vec & repeated_vector,
                               const char * link,
                               const arma::vec & beta_vector,
                               const arma::vec & mu_vector,
                               const arma::vec & eta_vector,
                               const arma::vec & alpha_vector,
                               const arma::vec & weights_vector);

arma::vec update_beta_robust_or(const arma::vec & y_vector,
                                const arma::mat & model_matrix,
                                const arma::vec & id_vector,
                                const arma::vec & repeated_vector,
                                const char * link,
                                const arma::vec & beta_vector,
                                const arma::vec & mu_vector,
                                const arma::vec & eta_vector,
                                const arma::vec & alpha_vector,
                                const arma::vec & weights_vector) ;

arma::vec update_beta_empirical_or(const arma::vec & y_vector,
                                   const arma::mat & model_matrix,
                                   const arma::vec & id_vector,
                                   const arma::vec & repeated_vector,
                                   const char * link,
                                   const arma::vec & beta_vector,
                                   const arma::vec & mu_vector,
                                   const arma::vec & eta_vector,
                                   const arma::vec & alpha_vector,
                                   const arma::vec & weights_vector);

arma::vec update_beta_jeffreys_or(const arma::vec & y_vector,
                                  const arma::mat & model_matrix,
                                  const arma::vec & id_vector,
                                  const arma::vec & repeated_vector,
                                  const char * link,
                                  const arma::vec & beta_vector,
                                  const arma::vec & mu_vector,
                                  const arma::vec & eta_vector,
                                  const arma::vec & alpha_vector,
                                  const arma::vec & weights_vector,
                                  const double & jeffreys_power);

arma::vec update_beta_or(const arma::vec & y_vector,
                         const arma::mat & model_matrix,
                         const arma::vec & id_vector,
                         const arma::vec & repeated_vector,
                         const char * link,
                         const arma::vec & beta_vector,
                         const arma::vec & mu_vector,
                         const arma::vec & eta_vector,
                         const arma::vec & alpha_vector,
                         const char * type,
                         const arma::vec & weights_vector,
                         const double & jeffreys_power);


#endif
