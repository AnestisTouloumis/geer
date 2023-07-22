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
                                      const arma::vec & alpha_vector) {
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
    arma::mat weight_matrix_i =
      get_weight_matrix_inverse_or(mu_vector_i,
                           get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                                            repeated_max,
                                                            alpha_vector));
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(d_matrix_i) * weight_matrix_i;
    arma::vec u_vector_i =
      t_d_matrix_weight_matrix_inverse_i * s_vector(id_vector_i);
    u_vector += u_vector_i;
    naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
    meat_matrix += u_vector_i * trans(u_vector_i);
  }
  arma::mat naive_matrix =
    arma::inv(naive_matrix_inverse, arma::inv_opts::allow_approx);
  arma::mat robust_matrix = naive_matrix * meat_matrix * naive_matrix;
  double obs_no_total = model_matrix.n_rows;
  double kappa_const = ((obs_no_total - 1)/(obs_no_total - params_no)) *
    (sample_size / (sample_size - 1));
  double delta_const = params_no / (sample_size - params_no);
  if (delta_const > 0.5) delta_const = 0.5;
  double phi_const = arma::trace(naive_matrix * meat_matrix)/ params_no;
  if (phi_const < 1.0) phi_const = 1.0;
  arma::mat bc_matrix =
    kappa_const * robust_matrix + delta_const * phi_const * naive_matrix;
  Rcpp::List ans;
  ans["naive_covariance"] = naive_matrix;
  ans["robust_covariance"] = robust_matrix;
  ans["bc_covariance"] = bc_matrix;
  return ans;
}
//==============================================================================

