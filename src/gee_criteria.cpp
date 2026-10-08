#include <stdexcept>
#include "link_functions.h"
#include "variance_functions.h"
#include "working_covariance_cc.h"
#include "cluster_utils.h"
#include "working_covariance_or.h"
#include "utils.h"

//============================ information matrix - independence =============
// [[Rcpp::export]]
arma::mat get_information_matrix_independence(const arma::mat& model_matrix,
                                              const arma::vec& id_vector,
                                              const char* link,
                                              const char* family,
                                              const arma::vec& mu_vector,
                                              const arma::vec& eta_vector,
                                              const double phi,
                                              const arma::vec& weights_vector) {
  const arma::uword params_no = model_matrix.n_cols;
  const arma::vec delta_vector = geer::link_derivative_1(link, eta_vector);
  const arma::vec variance_vector = geer::variance_function(family, mu_vector);
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::mat result(params_no, params_no, arma::fill::zeros);
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
      const arma::vec scale_i =
        weights_vector.subvec(first_row, last_row) / variance_vector.subvec(first_row, last_row);
      arma::mat d_scaled_i = d_matrix_i;
      d_scaled_i.each_col() %= scale_i;
      result += d_matrix_i.t() * d_scaled_i;
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_information_matrix_independence", cluster_index, id_vector[cluster.start], e);
    }
  }
  return result / phi;
}
//==============================================================================


//============================ get_sc_criteria =================================
// [[Rcpp::export]]
Rcpp::List get_working_covariance_criteria_cc(const arma::vec& y_vector,
                                              const arma::vec& id_vector,
                                              const arma::vec& repeated_vector,
                                              const char* family,
                                              const arma::vec& mu_vector,
                                              const char* correlation_structure,
                                              const arma::vec& alpha_vector,
                                              const double phi,
                                              const arma::vec& weights_vector) {
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  double sc_criterion = 0.0;
  double sum_log_det = 0.0;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      const arma::vec s_vector_i = s_vector.subvec(first_row, last_row);
      const arma::mat v_matrix_i =
        get_v_matrix_cc(family,
                        mu_vector.subvec(first_row, last_row),
                        repeated_vector.subvec(first_row, last_row),
                        phi,
                        correlation_matrix,
                        weights_vector.subvec(first_row, last_row));
      // One Cholesky factorization gives both s' V^{-1} s and log det(V).
      arma::mat chol_upper;
      if (!arma::chol(chol_upper, v_matrix_i)) {
        throw std::runtime_error(
          "working covariance V_i is not positive definite -- check correlation "
          "parameters and model specification.");
      }
      const arma::vec whitened_vector_i =
        arma::solve(arma::trimatl(chol_upper.t()), s_vector_i,
                    arma::solve_opts::no_approx);
      sc_criterion += arma::dot(whitened_vector_i, whitened_vector_i);
      sum_log_det += 2.0 * arma::accu(arma::log(chol_upper.diag()));
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_working_covariance_criteria_cc", cluster_index, id_vector[cluster.start], e);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("sc") = sc_criterion,
    Rcpp::Named("gp") = -0.5 * (sc_criterion + sum_log_det)
  );
}
//==============================================================================


//============================ sc criteria with odds ratios ====================
// [[Rcpp::export]]
Rcpp::List get_working_covariance_criteria_or(const arma::vec& y_vector,
                                              const arma::vec& id_vector,
                                              const arma::vec& repeated_vector,
                                              const arma::vec& mu_vector,
                                              const arma::vec& alpha_vector,
                                              const arma::vec& weights_vector) {
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  double sc_criterion = 0.0;
  double sum_log_det = 0.0;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      const arma::vec s_vector_i = s_vector.subvec(first_row, last_row);
      const arma::vec odds_ratios_vector_i =
        get_subject_specific_odds_ratios(repeated_vector.subvec(first_row, last_row),
                                         repeated_max,
                                         alpha_vector);
      const arma::mat v_matrix_i =
        get_v_matrix_or(mu_vector.subvec(first_row, last_row),
                        odds_ratios_vector_i,
                        weights_vector.subvec(first_row, last_row));
      // One Cholesky factorization gives both s' V^{-1} s and log det(V).
      arma::mat chol_upper;
      if (!arma::chol(chol_upper, v_matrix_i)) {
        throw std::runtime_error(
          "working covariance V_i is not positive definite -- check odds-ratio "
          "parameters and model specification.");
      }
      const arma::vec whitened_vector_i =
        arma::solve(arma::trimatl(chol_upper.t()), s_vector_i,
                    arma::solve_opts::no_approx);
      sc_criterion += arma::dot(whitened_vector_i, whitened_vector_i);
      sum_log_det += 2.0 * arma::accu(arma::log(chol_upper.diag()));
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_working_covariance_criteria_or", cluster_index, id_vector[cluster.start], e);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("sc") = sc_criterion,
    Rcpp::Named("gp") = -0.5 * (sc_criterion + sum_log_det)
  );
}
//==============================================================================
