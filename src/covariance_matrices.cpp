#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link_functions.h"
#include "nuisance_quantities_cc.h"
#include "nuisance_quantities_or.h"
#include "clusterutils.h"
using namespace Rcpp;


//============================ covariance matrices -- cc =======================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List get_covariance_matrices_cc(const arma::vec & y_vector,
                                      const arma::mat & model_matrix,
                                      const arma::vec & id_vector,
                                      const arma::vec & repeated_vector,
                                      const arma::vec & weights_vector,
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
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec s_vector = y_vector - mu_vector;
  arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure,
                           alpha_vector,
                           max(repeated_vector));
  auto clusters = clusters_from_sorted_id(id_vector);
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    arma::mat d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector.subvec(a, b);
    arma::mat v_matrix_i = get_v_matrix_cc(family,
                                           mu_vector.subvec(a, b),
                                           repeated_vector.subvec(a, b),
                                           phi,
                                           correlation_matrix,
                                           weights_vector.subvec(a, b));
    arma::mat d_matrix_trans_v_matrix_inverse_i =
      arma::solve(v_matrix_i, d_matrix_i, arma::solve_opts::likely_sympd).t();
    arma::vec u_vector_i =
      d_matrix_trans_v_matrix_inverse_i * s_vector.subvec(a, b);
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    meat_matrix += u_vector_i * u_vector_i.t();
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
  double ksi = arma::trace(naive_matrix_meat_matrix)/ params_no;
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


//============================ covariance matrices -- or =======================
// [[Rcpp::export()]]
Rcpp::List get_covariance_matrices_or(const arma::vec & y_vector,
                                      const arma::mat & model_matrix,
                                      const arma::vec & id_vector,
                                      const arma::vec & repeated_vector,
                                      const arma::vec & weights_vector,
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
  auto clusters = clusters_from_sorted_id(id_vector);
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    arma::mat d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector.subvec(a, b);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector.subvec(a, b),
                                       repeated_max,
                                       alpha_vector);
    arma::mat v_matrix_i = get_v_matrix_or(mu_vector.subvec(a, b),
                                           odds_ratios_vector_i,
                                           weights_vector.subvec(a, b));
    arma::mat d_matrix_trans_v_matrix_inverse_i =
      solve(v_matrix_i, d_matrix_i, arma::solve_opts::likely_sympd).t();
    arma::vec u_vector_i =
      d_matrix_trans_v_matrix_inverse_i * s_vector.subvec(a, b);
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    meat_matrix += u_vector_i * u_vector_i.t();
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
  double ksi = arma::trace(naive_matrix_meat_matrix)/ params_no;
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
