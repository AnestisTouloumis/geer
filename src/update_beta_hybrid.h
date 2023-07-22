#ifndef UPDATEBETAHYBRID_H
#define UPDATEBETAHYBRID_H

arma::vec update_beta_gee_hybrid(const arma::vec & y_vector,
                                 const arma::mat & model_matrix,
                                 const arma::vec & id_vector,
                                 const arma::vec & repeated_vector,
                                 const char * link,
                                 const arma::vec & beta_vector,
                                 const arma::vec & mu_vector,
                                 const arma::vec & eta_vector,
                                 const arma::vec & alpha_vector);

arma::vec update_beta_gee_hybrid(const arma::vec & y_vector,
                                 const arma::mat & model_matrix,
                                 const arma::vec & id_vector,
                                 const arma::vec & repeated_vector,
                                 const char * link,
                                 const arma::vec & beta_vector,
                                 const arma::vec & mu_vector,
                                 const arma::vec & eta_vector,
                                 const arma::vec & alpha_vector);

arma::vec update_beta_naive_hybrid(const arma::vec & y_vector,
                                   const arma::mat & model_matrix,
                                   const arma::vec & id_vector,
                                   const arma::vec & repeated_vector,
                                   const char * link,
                                   const arma::vec & beta_vector,
                                   const arma::vec & mu_vector,
                                   const arma::vec & eta_vector,
                                   const arma::vec & alpha_vector);

arma::vec update_beta_robust_hybrid(const arma::vec & y_vector,
                                    const arma::mat & model_matrix,
                                    const arma::vec & id_vector,
                                    const arma::vec & repeated_vector,
                                    const char * link,
                                    const arma::vec & beta_vector,
                                    const arma::vec & mu_vector,
                                    const arma::vec & eta_vector,
                                    const arma::vec & alpha_vector) ;

arma::vec update_beta_empirical_hybrid(const arma::vec & y_vector,
                                       const arma::mat & model_matrix,
                                       const arma::vec & id_vector,
                                       const arma::vec & repeated_vector,
                                       const char * link,
                                       const arma::vec & beta_vector,
                                       const arma::vec & mu_vector,
                                       const arma::vec & eta_vector,
                                       const arma::vec & alpha_vector);

arma::vec update_beta_jeffreys_hybrid(const arma::vec & y_vector,
                                      const arma::mat & model_matrix,
                                      const arma::vec & id_vector,
                                      const arma::vec & repeated_vector,
                                      const char * link,
                                      const arma::vec & beta_vector,
                                      const arma::vec & mu_vector,
                                      const arma::vec & eta_vector,
                                      const arma::vec & alpha_vector);

arma::vec update_beta_hybrid(const arma::vec & y_vector,
                             const arma::mat & model_matrix,
                             const arma::vec & id_vector,
                             const arma::vec & repeated_vector,
                             const char * link,
                             const arma::vec & beta_vector,
                             const arma::vec & mu_vector,
                             const arma::vec & eta_vector,
                             const arma::vec & alpha_vector,
                             const char * type);

#endif
