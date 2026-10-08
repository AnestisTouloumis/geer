#include "covariance_matrices.h"
#include "link_functions.h"
#include "working_covariance_cc.h"
#include "cluster_utils.h"
#include "working_covariance_or.h"
#include "utils.h"


// ========================== shared sandwich + BC correction ==================
static Rcpp::List compute_sandwich(arma::mat& information_matrix,
                                   const arma::mat& meat_matrix,
                                   const double sample_size,
                                   const double obs_no_total,
                                   const arma::uword params_no) {
  symmetrize_if_close(information_matrix, 1e-10);
  const arma::mat rhs =
    arma::join_rows(arma::eye(params_no, params_no), meat_matrix);
  const arma::mat solved = solve_chol_or_lu_mat(information_matrix, rhs);
  const arma::mat naive_covariance = solved.cols(0, params_no - 1);
  const arma::mat naive_covariance_meat_matrix = solved.cols(params_no,
                                                         2 * params_no - 1);
  arma::mat robust_matrix =
    solve_chol_or_lu_mat(information_matrix, naive_covariance_meat_matrix.t());
  symmetrize_if_close(robust_matrix, 1e-10);
  if (!(sample_size > static_cast<double>(params_no)) ||
      !(obs_no_total > static_cast<double>(params_no))) {
    arma::mat bc_undefined(params_no, params_no);
    bc_undefined.fill(NA_REAL);
    return Rcpp::List::create(
      Rcpp::Named("naive_covariance")  = naive_covariance,
      Rcpp::Named("robust_covariance") = robust_matrix,
      Rcpp::Named("bc_covariance")     = bc_undefined
    );
  }
  const double mbn_kappa =
    ((obs_no_total - 1.0) / (obs_no_total - static_cast<double>(params_no))) *
    (sample_size / (sample_size - 1.0));
  double mbn_lambda =
    static_cast<double>(params_no) /
      (sample_size - static_cast<double>(params_no));
  if (mbn_lambda > 0.5) mbn_lambda = 0.5;
  double mbn_xi =
    arma::trace(naive_covariance_meat_matrix) / static_cast<double>(params_no);
  if (mbn_xi < 1.0) mbn_xi = 1.0;
  arma::mat bc_matrix =
    mbn_kappa * robust_matrix + mbn_lambda * mbn_xi * naive_covariance;
  symmetrize_if_close(bc_matrix, 1e-10);
  return Rcpp::List::create(
    Rcpp::Named("naive_covariance")  = naive_covariance,
    Rcpp::Named("robust_covariance") = robust_matrix,
    Rcpp::Named("bc_covariance")     = bc_matrix
  );
}
//==============================================================================


//============================ covariance matrices -- cc =======================
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
                                      const double phi) {
  const auto clusters = clusters_from_sorted_id(id_vector);
  const double sample_size  = static_cast<double>(clusters.size());
  const arma::uword params_no = model_matrix.n_cols;
  const double obs_no_total = static_cast<double>(model_matrix.n_rows);
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = geer::link_derivative_1(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
  arma::mat d_matrix_i;
  arma::mat v_matrix_i;
  arma::mat v_matrix_inverse_d_matrix_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::vec u_vector_i;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      d_matrix_i = model_matrix.rows(first_row, last_row);
      d_matrix_i.each_col() %= delta_vector.subvec(first_row, last_row);
      v_matrix_i = get_v_matrix_cc(family,
                                   mu_vector.subvec(first_row, last_row),
                                   repeated_vector.subvec(first_row, last_row),
                                   phi,
                                   correlation_matrix,
                                   weights_vector.subvec(first_row, last_row));
      v_matrix_inverse_d_matrix_i =
        solve_chol_or_lu_mat(v_matrix_i, d_matrix_i);
      d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
      u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector.subvec(first_row, last_row);
      information_matrix += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      meat_matrix += u_vector_i * u_vector_i.t();
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_covariance_matrices_cc", cluster_index, id_vector[cluster.start], e);
    }
  }
  return compute_sandwich(information_matrix, meat_matrix,
                          sample_size, obs_no_total, params_no);
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
  const auto clusters = clusters_from_sorted_id(id_vector);
  const double sample_size  = static_cast<double>(clusters.size());
  const arma::uword params_no = model_matrix.n_cols;
  const double obs_no_total = static_cast<double>(model_matrix.n_rows);
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = geer::link_derivative_1(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  arma::mat d_matrix_i;
  arma::mat v_matrix_i;
  arma::mat v_matrix_inverse_d_matrix_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::vec u_vector_i;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      d_matrix_i = model_matrix.rows(first_row, last_row);
      d_matrix_i.each_col() %= delta_vector.subvec(first_row, last_row);
      const arma::vec odds_ratios_vector_i =
        get_subject_specific_odds_ratios(repeated_vector.subvec(first_row, last_row),
                                         repeated_max,
                                         alpha_vector);
      v_matrix_i = get_v_matrix_or(mu_vector.subvec(first_row, last_row),
                                   odds_ratios_vector_i,
                                   weights_vector.subvec(first_row, last_row));
      v_matrix_inverse_d_matrix_i =
        solve_chol_or_lu_mat(v_matrix_i, d_matrix_i);
      d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
      u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector.subvec(first_row, last_row);
      information_matrix += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      meat_matrix += u_vector_i * u_vector_i.t();
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_covariance_matrices_or", cluster_index, id_vector[cluster.start], e);
    }
  }
  return compute_sandwich(information_matrix, meat_matrix,
                          sample_size, obs_no_total, params_no);
}
//==============================================================================
