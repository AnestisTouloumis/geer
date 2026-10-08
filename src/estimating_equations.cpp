#include <cmath>
#include "link_functions.h"
#include "working_covariance_cc.h"
#include "cluster_utils.h"
#include "working_covariance_or.h"
#include "utils.h"

namespace {
inline void validate_estimating_equations_inputs(const arma::vec& y_vector,
                                                 const arma::mat& model_matrix,
                                                 const arma::vec& id_vector,
                                                 const arma::vec& repeated_vector,
                                                 const arma::vec& weights_vector,
                                                 const arma::vec& mu_vector,
                                                 const arma::vec& eta_vector) {
  const arma::uword n = y_vector.n_elem;
  if (n == 0) {
    Rcpp::stop("Input vectors must not be empty.");
  }
  if (model_matrix.n_rows != n) {
    Rcpp::stop("'model_matrix' must have the same number of rows as 'y_vector'.");
  }
  if (id_vector.n_elem != n || repeated_vector.n_elem != n ||
      weights_vector.n_elem != n || mu_vector.n_elem != n ||
      eta_vector.n_elem != n) {
    Rcpp::stop("All observation-level inputs must have the same length.");
  }
}
}

//============================ estimating equations - cc =======================
// [[Rcpp::export]]
arma::vec estimating_equations_gee_cc(const arma::vec& y_vector,
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
  validate_estimating_equations_inputs(y_vector,
                             model_matrix,
                             id_vector,
                             repeated_vector,
                             weights_vector,
                             mu_vector,
                             eta_vector);

  if (!R_FINITE(phi) || phi <= 0.0) {
    Rcpp::stop("'phi' must be finite and positive.");
  }
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
  const arma::vec delta_vector = geer::link_derivative_1(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::vec result(params_no, arma::fill::zeros);
  arma::mat d_matrix_i;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;

      d_matrix_i = model_matrix.rows(first_row, last_row);
      d_matrix_i.each_col() %= delta_vector.subvec(first_row, last_row);

      const arma::mat v_matrix_i =
        get_v_matrix_cc(family,
                        mu_vector.subvec(first_row, last_row),
                        repeated_vector.subvec(first_row, last_row),
                        phi,
                        correlation_matrix,
                        weights_vector.subvec(first_row, last_row));

      result += d_matrix_i.t() *
        solve_chol_or_lu_vec(v_matrix_i, s_vector.subvec(first_row, last_row));
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("estimating_equations_gee_cc", cluster_index, id_vector[cluster.start], e);
    }
  }
  return result;
}
//==============================================================================


//============================ estimating equations - or =======================
// [[Rcpp::export]]
arma::vec estimating_equations_gee_or(const arma::vec& y_vector,
                                      const arma::mat& model_matrix,
                                      const arma::vec& id_vector,
                                      const arma::vec& repeated_vector,
                                      const arma::vec& weights_vector,
                                      const char* link,
                                      const arma::vec& mu_vector,
                                      const arma::vec& eta_vector,
                                      const arma::vec& alpha_vector) {
  validate_estimating_equations_inputs(y_vector,
                             model_matrix,
                             id_vector,
                             repeated_vector,
                             weights_vector,
                             mu_vector,
                             eta_vector);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec delta_vector = geer::link_derivative_1(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::vec result(params_no, arma::fill::zeros);
  arma::mat d_matrix_i;
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
      const arma::mat v_matrix_i =
        get_v_matrix_or(mu_vector.subvec(first_row, last_row),
                        odds_ratios_vector_i,
                        weights_vector.subvec(first_row, last_row));

      result += d_matrix_i.t() *
        solve_chol_or_lu_vec(v_matrix_i, s_vector.subvec(first_row, last_row));
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("estimating_equations_gee_or", cluster_index, id_vector[cluster.start], e);
    }
  }
  return result;
}
//==============================================================================
