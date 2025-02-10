#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link_functions.h"
#include "utils.h"
#include "nuisance_quantities_oddsratio.h"
using namespace Rcpp;


//============================ covariance matrices with odds ratios ============
// [[Rcpp::export()]]
Rcpp::List get_covariance_matrices_or(const arma::vec & y_vector,
                                      const arma::mat & model_matrix,
                                      const arma::vec & id_vector,
                                      const arma::vec & repeated_vector,
                                      const char * link,
                                      const arma::vec & mu_vector,
                                      const arma::vec & eta_vector,
                                      const arma::vec & alpha_vector,
                                      const arma::vec & weights_vector) {
  double sample_size = max(id_vector);
  double params_no = model_matrix.n_cols;
  int repeated_max = max(repeated_vector);
  arma::vec u_vector = arma::zeros(params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::mat meat_matrix = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec s_vector = y_vector - mu_vector;
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec mu_vector_i = mu_vector(id_vector_i);
    arma::mat d_matrix_i =
      arma::diagmat(delta_vector(id_vector_i)) *
      model_matrix.rows(id_vector_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max,
                                       alpha_vector);
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(
        solve(get_weight_matrix_or(mu_vector(id_vector_i),
                                   odds_ratios_vector_i,
                                   weights_vector(id_vector_i)),
                                   d_matrix_i));
    arma::vec u_vector_i =
      t_d_matrix_weight_matrix_inverse_i * s_vector(id_vector_i);
    naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
    meat_matrix += u_vector_i * trans(u_vector_i);
  }
  arma::mat naive_matrix = arma::pinv(naive_matrix_inverse);
  arma::mat naive_matrix_meat_matrix = solve(naive_matrix_inverse, meat_matrix);
  arma::mat robust_matrix = solve(naive_matrix_inverse,
                                  trans(naive_matrix_meat_matrix));
  double obs_no_total = model_matrix.n_rows;
  double kappa = ((obs_no_total - 1)/(obs_no_total - params_no)) *
    (sample_size / (sample_size - 1));
  double lambda = params_no / (sample_size - params_no);
  if (lambda > 0.5) lambda = 0.5;
  double ksi = arma::trace(naive_matrix * meat_matrix)/ params_no;
  if (ksi < 1.0) ksi = 1.0;
  arma::mat bc_matrix =
    kappa * robust_matrix + lambda * ksi * naive_matrix;
  Rcpp::List ans;
  ans["naive_covariance"] = naive_matrix;
  ans["robust_covariance"] = robust_matrix;
  ans["bc_covariance"] = bc_matrix;
  return ans;
}
//==============================================================================

