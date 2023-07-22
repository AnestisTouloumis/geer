#ifndef NUISANCEQUANTITIESODDSRATIO_H
#define NUISANCEQUANTITIESODDSRATIO_H

Rcpp::NumericVector get_marginalized_odds_ratios(const arma::vec & response_vector,
                                                 const arma::vec & id_vector,
                                                 const arma::vec & repeated_vector,
                                                 const double &  adding_constant,
                                                 Rcpp::String or_structure);

arma::vec get_subject_specific_odds_ratios(arma::vec repeated_vector_i,
                                           int cluster_size_max,
                                           arma::vec odds_ratios_vector);

double get_bivariate_distribution(const double & row_prob,
                                  const double & col_prob,
                                  const double & odds_ratio);

arma::mat get_weight_matrix_or(const arma::vec & mu_vector,
                               const arma::vec & odds_ratios_vector);

arma::mat get_weight_matrix_inverse_or(const arma::vec & mu_vector,
                                       const arma::vec & odds_ratios_vector);

double get_bivariate_distribution_murow(const double & row_prob,
                                        const double & col_prob,
                                        const double & odds_ratio);

double get_bivariate_distribution_murow2(const double & row_prob,
                                         const double & col_prob,
                                         const double & odds_ratio);

double get_bivariate_distribution_murowcol(double row_prob,
                                           double col_prob,
                                           double odds_ratio);

arma::mat get_g_matrix(const arma::vec & mu_vector,
                       const arma::vec & odds_ratios_vector);

arma::mat get_bivariate_distribution_mu(const arma::vec & mu_vector,
                                        const arma::vec & odds_ratios_vector);

arma::mat get_g_matrix_mu(const arma::vec & mu_vector,
                          const arma::vec & odds_ratios_vector);

arma::mat get_weight_matrix_mu_or(const arma::vec & mu_vector,
                                  const arma::vec & odds_ratios_vector);

#endif
