#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link_functions.h"
#include "variance_functions.h"
#include "nuisance_quantities.h"
#include "utils.h"
using namespace Rcpp;


//============================ covariance matrices =============================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
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
                                   const double & phi) {
  double sample_size = max(id_vector);
  double params_no = model_matrix.n_cols;
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::mat meat_matrix = arma::zeros(params_no, params_no);
  arma::vec delta_star_vector = mueta(link, eta_vector);
  arma::vec s_vector = y_vector - mu_vector;
  arma::mat correlation_matrix_inverse_full =
    get_correlation_matrix_inverse(correlation_structure,
                                   alpha_vector,
                                   max(repeated_vector));
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::mat d_matrix_i =
      arma::diagmat(delta_star_vector(id_vector_i)) *
      model_matrix.rows(id_vector_i);
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(d_matrix_i) *
      get_weight_matrix(family, mu_vector(id_vector_i), repeated_vector(id_vector_i),
                              phi, correlation_matrix_inverse_full);
    arma::vec u_vector_i =
      t_d_matrix_weight_matrix_inverse_i * s_vector(id_vector_i);
    naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
    meat_matrix += u_vector_i * trans(u_vector_i);
  }
  arma::mat naive_matrix =
    arma::inv(naive_matrix_inverse, arma::inv_opts::allow_approx);
  arma::mat robust_matrix = naive_matrix * meat_matrix * naive_matrix;
  double obs_no_total = model_matrix.n_rows;
  double kappa_const =
    ((obs_no_total - 1)/(obs_no_total - params_no)) * (sample_size / (sample_size - 1));
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
