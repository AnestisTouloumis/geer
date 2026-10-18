#include "small_sample_covariances.h"
#include "cluster_utils.h"
#include "link_functions.h"
#include "utils.h"
#include "working_covariance_cc.h"
#include "working_covariance_or.h"


namespace {
//============================ per-cluster quantities ==========================
struct ClusterQuantities {
  double id_value;
  arma::mat d_matrix;
  arma::mat v_inverse_d_matrix;
  arma::vec residuals;
  arma::vec sd_vector;
};
//==============================================================================


//============================ small-sample covariance estimators ==============
Rcpp::List get_small_sample_covariances(
  const std::vector<ClusterQuantities>& clusters,
  const arma::uword params_no,
  const arma::uword repeated_max) {
  const arma::uword clusters_no = clusters.size();
  const double sample_size = static_cast<double>(clusters_no);
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  bool complete_clusters = true;
  for (const auto& cluster : clusters) {
    const arma::vec u_vector = cluster.v_inverse_d_matrix.t() * cluster.residuals;
    information_matrix += cluster.v_inverse_d_matrix.t() * cluster.d_matrix;
    meat_matrix += u_vector * u_vector.t();
    if (cluster.residuals.n_elem != repeated_max) complete_clusters = false;
  }
  symmetrize_if_close(information_matrix, 1e-10);
  const arma::mat naive_covariance = solve_chol_or_lu_mat(
    information_matrix, arma::eye(params_no, params_no),
    "get_small_sample_covariances");

  arma::mat meat_kc(params_no, params_no, arma::fill::zeros);
  arma::mat meat_md(params_no, params_no, arma::fill::zeros);
  arma::mat meat_fg(params_no, params_no, arma::fill::zeros);
  arma::mat meat_fz(params_no, params_no, arma::fill::zeros);
  arma::mat pooled_pan(repeated_max, repeated_max, arma::fill::zeros);
  arma::mat pooled_wl(repeated_max, repeated_max, arma::fill::zeros);
  arma::mat pooled_wb(repeated_max, repeated_max, arma::fill::zeros);
  for (arma::uword cluster_index = 0; cluster_index < clusters_no; ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::mat& d_matrix = cluster.d_matrix;
      const arma::mat& v_inverse_d_matrix = cluster.v_inverse_d_matrix;
      const arma::vec& residuals = cluster.residuals;
      const arma::uword cluster_size = residuals.n_elem;
      const arma::vec u_vector = v_inverse_d_matrix.t() * residuals;
      const arma::mat h_matrix = d_matrix * naive_covariance * v_inverse_d_matrix.t();
      const arma::mat md_matrix =
        arma::inv(arma::eye(cluster_size, cluster_size) - h_matrix);
      const arma::mat kc_matrix = arma::real(arma::sqrtmat(md_matrix));
      const arma::vec md_residuals = md_matrix * residuals;
      const arma::vec kc_residuals = kc_matrix * residuals;
      const arma::vec md_u_vector = v_inverse_d_matrix.t() * md_residuals;
      const arma::vec kc_u_vector = v_inverse_d_matrix.t() * kc_residuals;
      meat_md += md_u_vector * md_u_vector.t();
      meat_kc += kc_u_vector * kc_u_vector.t();

      const arma::vec leverage_vector =
        arma::diagvec(v_inverse_d_matrix.t() * d_matrix * naive_covariance);
      const arma::vec fg_vector =
        1.0 / arma::sqrt(1.0 - arma::clamp(leverage_vector, -arma::datum::inf, 0.75));
      const arma::vec fg_u_vector = fg_vector % u_vector;
      meat_fg += fg_u_vector * fg_u_vector.t();

      const arma::mat inner_matrix =
        residuals * residuals.t() -
        d_matrix * naive_covariance * (meat_matrix - u_vector * u_vector.t()) *
          naive_covariance * d_matrix.t();
      const arma::mat fz_matrix = v_inverse_d_matrix.t() * md_matrix;
      meat_fz += fz_matrix * inner_matrix * fz_matrix.t();

      if (complete_clusters) {
        const arma::vec pan_residuals = residuals / cluster.sd_vector;
        const arma::vec wl_residuals = md_residuals / cluster.sd_vector;
        const arma::vec wb_residuals = kc_residuals / cluster.sd_vector;
        pooled_pan += pan_residuals * pan_residuals.t();
        pooled_wl += wl_residuals * wl_residuals.t();
        pooled_wb += wb_residuals * wb_residuals.t();
      }
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_small_sample_covariances", cluster_index,
                                   cluster.id_value, e);
    }
  }

  arma::mat na_matrix(params_no, params_no);
  na_matrix.fill(NA_REAL);
  arma::mat meat_pan = na_matrix;
  arma::mat meat_wl = na_matrix;
  arma::mat meat_wb = na_matrix;
  if (complete_clusters) {
    meat_pan.zeros();
    meat_wl.zeros();
    meat_wb.zeros();
    pooled_pan /= sample_size;
    pooled_wl /= sample_size;
    pooled_wb /= sample_size;
    for (const auto& cluster : clusters) {
      arma::mat scaled_matrix = cluster.v_inverse_d_matrix;
      scaled_matrix.each_col() %= cluster.sd_vector;
      meat_pan += scaled_matrix.t() * pooled_pan * scaled_matrix;
      meat_wl += scaled_matrix.t() * pooled_wl * scaled_matrix;
      meat_wb += scaled_matrix.t() * pooled_wb * scaled_matrix;
    }
  }

  auto sandwich = [&naive_covariance](const arma::mat& meat) {
    arma::mat covariance = naive_covariance * meat * naive_covariance;
    symmetrize_if_close(covariance, 1e-10);
    return covariance;
  };
  const arma::mat kc_covariance = sandwich(meat_kc);
  const arma::mat md_covariance = sandwich(meat_md);
  const arma::mat pan_covariance = sandwich(meat_pan);
  arma::mat gst_covariance = na_matrix;
  if (complete_clusters && sample_size > static_cast<double>(params_no)) {
    gst_covariance =
      sample_size / (sample_size - static_cast<double>(params_no)) * pan_covariance;
  }
  return Rcpp::List::create(
    Rcpp::Named("kauermann-carroll")   = kc_covariance,
    Rcpp::Named("mancl-derouen")       = md_covariance,
    Rcpp::Named("ford-westgate")       = 0.5 * (kc_covariance + md_covariance),
    Rcpp::Named("fay-graubard")        = sandwich(meat_fg),
    Rcpp::Named("fan-zhang-zhang")     = sandwich(meat_fz),
    Rcpp::Named("pan")                 = pan_covariance,
    Rcpp::Named("gosho-sato-takeuchi") = gst_covariance,
    Rcpp::Named("wang-long")           = sandwich(meat_wl),
    Rcpp::Named("westgate-burchett")   = sandwich(meat_wb)
  );
}
//==============================================================================
}  // namespace


//============================ small-sample covariances -- cc ==================
// [[Rcpp::export]]
Rcpp::List get_small_sample_covariances_cc(const arma::vec& y_vector,
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
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec delta_vector = link_derivative_1(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
  std::vector<ClusterQuantities> quantities;
  quantities.reserve(clusters.size());
  for (arma::uword cluster_index = 0; cluster_index < clusters.size(); ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      arma::mat d_matrix_i = model_matrix.rows(first_row, last_row);
      d_matrix_i.each_col() %= delta_vector.subvec(first_row, last_row);
      const arma::mat v_matrix_i =
        get_v_matrix_cc(family,
                        mu_vector.subvec(first_row, last_row),
                        repeated_vector.subvec(first_row, last_row),
                        phi,
                        correlation_matrix,
                        weights_vector.subvec(first_row, last_row));
      quantities.push_back({id_vector[first_row],
                            d_matrix_i,
                            solve_chol_or_lu_mat(v_matrix_i, d_matrix_i),
                            s_vector.subvec(first_row, last_row),
                            arma::sqrt(arma::diagvec(v_matrix_i))});
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_small_sample_covariances_cc", cluster_index,
                                   id_vector[cluster.start], e);
    }
  }
  return get_small_sample_covariances(quantities, model_matrix.n_cols, repeated_max);
}
//==============================================================================


//============================ small-sample covariances -- or ==================
// [[Rcpp::export]]
Rcpp::List get_small_sample_covariances_or(const arma::vec& y_vector,
                                           const arma::mat& model_matrix,
                                           const arma::vec& id_vector,
                                           const arma::vec& repeated_vector,
                                           const arma::vec& weights_vector,
                                           const char* link,
                                           const arma::vec& mu_vector,
                                           const arma::vec& eta_vector,
                                           const arma::vec& alpha_vector) {
  const auto clusters = clusters_from_sorted_id(id_vector);
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec delta_vector = link_derivative_1(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  std::vector<ClusterQuantities> quantities;
  quantities.reserve(clusters.size());
  for (arma::uword cluster_index = 0; cluster_index < clusters.size(); ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      arma::mat d_matrix_i = model_matrix.rows(first_row, last_row);
      d_matrix_i.each_col() %= delta_vector.subvec(first_row, last_row);
      const arma::vec odds_ratios_vector_i =
        get_subject_specific_odds_ratios(repeated_vector.subvec(first_row, last_row),
                                         repeated_max,
                                         alpha_vector);
      const arma::mat v_matrix_i =
        get_v_matrix_or(mu_vector.subvec(first_row, last_row),
                        odds_ratios_vector_i,
                        weights_vector.subvec(first_row, last_row));
      quantities.push_back({id_vector[first_row],
                            d_matrix_i,
                            solve_chol_or_lu_mat(v_matrix_i, d_matrix_i),
                            s_vector.subvec(first_row, last_row),
                            arma::sqrt(arma::diagvec(v_matrix_i))});
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_small_sample_covariances_or", cluster_index,
                                   id_vector[cluster.start], e);
    }
  }
  return get_small_sample_covariances(quantities, model_matrix.n_cols, repeated_max);
}
//==============================================================================
