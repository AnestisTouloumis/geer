#ifndef UPDATEBETA_H
#define UPDATEBETA_H

arma::vec update_beta_gee(const arma::vec & y_vector,
                          const arma::mat & model_matrix,
                          const arma::vec & id_vector,
                          const arma::vec & repeated_vector,
                          const char* link,
                          const char* family,
                          const arma::vec & beta_vector,
                          const arma::vec & mu_vector,
                          const arma::vec & eta_vector,
                          const char * correlation_structure,
                          const arma::vec & alpha_vector,
                          const double & phi);

arma::vec update_beta_naive(const arma::vec & y_vector,
                            const arma::mat & model_matrix,
                            const arma::vec & id_vector,
                            const arma::vec & repeated_vector,
                            const char* link,
                            const char* family,
                            const arma::vec & beta_vector,
                            const arma::vec & mu_vector,
                            const arma::vec & eta_vector,
                            const char * correlation_structure,
                            const arma::vec & alpha_vector,
                            const double & phi);

arma::vec update_beta_robust(const arma::vec & y_vector,
                             const arma::mat & model_matrix,
                             const arma::vec & id_vector,
                             const arma::vec & repeated_vector,
                             const char* link,
                             const char* family,
                             const arma::vec & beta_vector,
                             const arma::vec & mu_vector,
                             const arma::vec & eta_vector,
                             const char * correlation_structure,
                             const arma::vec & alpha_vector,
                             const double & phi) ;

arma::vec update_beta_empirical(const arma::vec & y_vector,
                                const arma::mat & model_matrix,
                                const arma::vec & id_vector,
                                const arma::vec & repeated_vector,
                                const char* link,
                                const char* family,
                                const arma::vec & beta_vector,
                                const arma::vec & mu_vector,
                                const arma::vec & eta_vector,
                                const char * correlation_structure,
                                const arma::vec & alpha_vector,
                                const double & phi);

arma::vec update_beta_jeffreys(const arma::vec & y_vector,
                               const arma::mat & model_matrix,
                               const arma::vec & id_vector,
                               const arma::vec & repeated_vector,
                               const char* link,
                               const char* family,
                               const arma::vec & beta_vector,
                               const arma::vec & mu_vector,
                               const arma::vec & eta_vector,
                               const char * correlation_structure,
                               const arma::vec & alpha_vector,
                               const double & phi);


arma::vec update_beta(const arma::vec & y_vector,
                      const arma::mat & model_matrix,
                      const arma::vec & id_vector,
                      const arma::vec & repeated_vector,
                      const char* link,
                      const char* family,
                      const arma::vec & beta_vector,
                      const arma::vec & mu_vector,
                      const arma::vec & eta_vector,
                      const char * correlation_structure,
                      const arma::vec & alpha_vector,
                      const double & phi,
                      const char* type);

#endif
