#include "cluster_utils.h"
#include "link_functions.h"
#include "utils.h"
#include "working_covariance_cc.h"
#include "working_covariance_or.h"
#include <cmath>


//============================ estimating equations - cc =======================
// [[Rcpp::export]]
arma::vec get_estimating_equations_cc(const arma::vec& y_vector,
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
  reset_lu_fallback_warning();
  check_cluster_inputs(id_vector,
                       {y_vector, repeated_vector, weights_vector, mu_vector,
                        eta_vector},
                       model_matrix);

  if (!R_FINITE(phi) || phi <= 0.0) {
    signal_input_error("'phi' must be a positive finite number");
  }
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::mat correlation_matrix =
    get_correlation_matrix_cc(correlation_structure, alpha_vector, repeated_max);
  const arma::vec delta_vector = link_derivative_1(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::vec result(params_no, arma::fill::zeros);
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const ClusterWorkspace workspace = build_cluster_data_cc(
        cluster, model_matrix, mu_vector, repeated_vector, weights_vector,
        delta_vector, s_vector, parse_family(family), phi, correlation_matrix);
      const ClusterData& data = workspace.data;
      const arma::mat& d_matrix_i = data.d_matrix;
      const arma::mat& v_matrix_i = data.v_matrix;
      result += d_matrix_i.t() *
                solve_chol_or_lu_vec(v_matrix_i, data.residuals);
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_estimating_equations_cc", cluster_index,
                                   id_vector[cluster.first_row], e);
    }
  }
  return result;
}
//==============================================================================


//============================ estimating equations - or =======================
// [[Rcpp::export]]
arma::vec get_estimating_equations_or(const arma::vec& y_vector,
                                      const arma::mat& model_matrix,
                                      const arma::vec& id_vector,
                                      const arma::vec& repeated_vector,
                                      const arma::vec& weights_vector,
                                      const char* link,
                                      const arma::vec& mu_vector,
                                      const arma::vec& eta_vector,
                                      const arma::vec& alpha_vector) {
  reset_lu_fallback_warning();
  check_cluster_inputs(id_vector,
                       {y_vector, repeated_vector, weights_vector, mu_vector,
                        eta_vector},
                       model_matrix);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec delta_vector = link_derivative_1(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::vec result(params_no, arma::fill::zeros);
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const ClusterWorkspace workspace = build_cluster_data_or(
        cluster, model_matrix, mu_vector, repeated_vector, weights_vector,
        delta_vector, s_vector, repeated_max, alpha_vector);
      const ClusterData& data = workspace.data;
      const arma::mat& d_matrix_i = data.d_matrix;
      const arma::mat& v_matrix_i = data.v_matrix;
      result += d_matrix_i.t() *
                solve_chol_or_lu_vec(v_matrix_i, data.residuals);
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_estimating_equations_or", cluster_index,
                                   id_vector[cluster.first_row], e);
    }
  }
  return result;
}
//==============================================================================
