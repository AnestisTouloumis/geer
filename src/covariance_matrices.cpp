#define ARMA_WARN_LEVEL 1
#include "covariance_matrices.h"
#include "link_functions.h"
#include "nuisance_quantities_cc.h"
#include "nuisance_quantities_or.h"
#include "clusterutils.h"
#include "utils.h"


//============================ covariance matrices -- cc =======================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List get_covariance_matrices_cc(const arma::vec& y_vector,
                                      const arma::mat& model_matrix,
                                      const arma::vec& id_vector,
                                      const arma::vec& repeated_vector,
                                      const arma::vec& weights_vector,
                                      const char* link,
                                      const char* family,
                                      const arma::vec& mu_vector,
                                      const arma::vec& eta_vector,
                                      const char* correlation_structure,
                                      const arma::vec& alpha_vector,
                                      const double& phi) {
  const double sample_size = arma::max(id_vector);
  const arma::uword params_no = model_matrix.n_cols;
  const double obs_no_total = static_cast<double>(model_matrix.n_rows);
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::mat naive_matrix_inverse(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = mueta(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure,
                           alpha_vector,
                           repeated_max);
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::mat d_matrix_i;
  arma::mat v_matrix_i;
  arma::mat v_matrix_inverse_d_matrix_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::vec u_vector_i(params_no, arma::fill::zeros);
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::uword m = cl.end - cl.start;
    if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
      d_matrix_i.set_size(m, params_no);
    }
    d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector.subvec(a, b);
    v_matrix_i = get_v_matrix_cc(family,
                                 mu_vector.subvec(a, b),
                                 repeated_vector.subvec(a, b),
                                 phi,
                                 correlation_matrix,
                                 weights_vector.subvec(a, b));
    if (v_matrix_inverse_d_matrix_i.n_rows != m ||
        v_matrix_inverse_d_matrix_i.n_cols != params_no) {
      v_matrix_inverse_d_matrix_i.set_size(m, params_no);
    }
    const bool ok_cluster = arma::solve(v_matrix_inverse_d_matrix_i,
                                        v_matrix_i,
                                        d_matrix_i,
                                        arma::solve_opts::likely_sympd);
    if (!ok_cluster) {
      Rcpp::stop("get_covariance_matrices_cc: failed to solve V_i^{-1} D_i.");
    }
    if (d_matrix_trans_v_matrix_inverse_i.n_rows != params_no ||
        d_matrix_trans_v_matrix_inverse_i.n_cols != m) {
      d_matrix_trans_v_matrix_inverse_i.set_size(params_no, m);
    }
    d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
    u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector.subvec(a, b);
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    meat_matrix += u_vector_i * u_vector_i.t();
  }
  symmetrize_if_close(naive_matrix_inverse, 1e-10);
  const arma::mat naive_matrix =
    solve_chol_or_lu_mat(naive_matrix_inverse, arma::eye(params_no, params_no));
  const arma::mat naive_matrix_meat_matrix =
    solve_chol_or_lu_mat(naive_matrix_inverse, meat_matrix);
  const arma::mat robust_matrix =
    solve_chol_or_lu_mat(naive_matrix_inverse, naive_matrix_meat_matrix.t());
  const double kappa =
    ((obs_no_total - 1.0) / (obs_no_total - static_cast<double>(params_no))) *
    (sample_size / (sample_size - 1.0));
  double lambda =
    static_cast<double>(params_no) / (sample_size - static_cast<double>(params_no));
  if (lambda > 0.5) {
    lambda = 0.5;
  }
  double ksi = arma::trace(naive_matrix_meat_matrix) / static_cast<double>(params_no);
  if (ksi < 1.0) {
    ksi = 1.0;
  }
  const arma::mat bc_matrix =
    kappa * robust_matrix + lambda * ksi * naive_matrix;
  return Rcpp::List::create(
    Rcpp::Named("naive_covariance") = naive_matrix,
    Rcpp::Named("robust_covariance") = robust_matrix,
    Rcpp::Named("bc_covariance") = bc_matrix
  );
}
//==============================================================================


//============================ covariance matrices -- or =======================
// [[Rcpp::export]]
Rcpp::List get_covariance_matrices_or(const arma::vec& y_vector,
                                      const arma::mat& model_matrix,
                                      const arma::vec& id_vector,
                                      const arma::vec& repeated_vector,
                                      const arma::vec& weights_vector,
                                      const char* link,
                                      const arma::vec& mu_vector,
                                      const arma::vec& eta_vector,
                                      const arma::vec& alpha_vector) {
  const double sample_size = arma::max(id_vector);
  const arma::uword params_no = model_matrix.n_cols;
  const double obs_no_total = static_cast<double>(model_matrix.n_rows);
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::mat naive_matrix_inverse(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = mueta(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::mat d_matrix_i;
  arma::mat v_matrix_i;
  arma::mat v_matrix_inverse_d_matrix_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::vec u_vector_i(params_no, arma::fill::zeros);
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::uword m = cl.end - cl.start;
    if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
      d_matrix_i.set_size(m, params_no);
    }
    d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector.subvec(a, b);
    const arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector.subvec(a, b),
                                       repeated_max,
                                       alpha_vector);
    v_matrix_i = get_v_matrix_or(mu_vector.subvec(a, b),
                                 odds_ratios_vector_i,
                                 weights_vector.subvec(a, b));
    if (v_matrix_inverse_d_matrix_i.n_rows != m ||
        v_matrix_inverse_d_matrix_i.n_cols != params_no) {
      v_matrix_inverse_d_matrix_i.set_size(m, params_no);
    }
    const bool ok_cluster = arma::solve(v_matrix_inverse_d_matrix_i,
                                        v_matrix_i,
                                        d_matrix_i,
                                        arma::solve_opts::likely_sympd);
    if (!ok_cluster) {
      Rcpp::stop("get_covariance_matrices_or: failed to solve V_i^{-1} D_i.");
    }
    if (d_matrix_trans_v_matrix_inverse_i.n_rows != params_no ||
        d_matrix_trans_v_matrix_inverse_i.n_cols != m) {
      d_matrix_trans_v_matrix_inverse_i.set_size(params_no, m);
    }
    d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
    u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector.subvec(a, b);
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    meat_matrix += u_vector_i * u_vector_i.t();
  }
  symmetrize_if_close(naive_matrix_inverse, 1e-10);
  const arma::mat naive_matrix =
    solve_chol_or_lu_mat(naive_matrix_inverse, arma::eye(params_no, params_no));
  const arma::mat naive_matrix_meat_matrix =
    solve_chol_or_lu_mat(naive_matrix_inverse, meat_matrix);
  const arma::mat robust_matrix =
    solve_chol_or_lu_mat(naive_matrix_inverse, naive_matrix_meat_matrix.t());
  const double kappa =
    ((obs_no_total - 1.0) / (obs_no_total - static_cast<double>(params_no))) *
    (sample_size / (sample_size - 1.0));
  double lambda =
    static_cast<double>(params_no) / (sample_size - static_cast<double>(params_no));
  if (lambda > 0.5) {
    lambda = 0.5;
  }
  double ksi = arma::trace(naive_matrix_meat_matrix) / static_cast<double>(params_no);
  if (ksi < 1.0) {
    ksi = 1.0;
  }
  const arma::mat bc_matrix =
    kappa * robust_matrix + lambda * ksi * naive_matrix;
  return Rcpp::List::create(
    Rcpp::Named("naive_covariance") = naive_matrix,
    Rcpp::Named("robust_covariance") = robust_matrix,
    Rcpp::Named("bc_covariance") = bc_matrix
  );
}
//==============================================================================
