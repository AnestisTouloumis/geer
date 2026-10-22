#include "small_sample_covariances.h"
#include "cluster_utils.h"
#include "link_functions.h"
#include "utils.h"
#include "working_covariance_cc.h"
#include "working_covariance_or.h"
#include <algorithm>
#include <stdexcept>
#include <vector>


namespace {
//============================ small-sample covariance estimators ==============
Rcpp::List get_small_sample_covariances(
  const std::vector<ClusterData>& clusters,
  const arma::uword params_no,
  const arma::uword repeated_max) {
  const arma::uword clusters_no = clusters.size();
  const double sample_size = static_cast<double>(clusters_no);
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  bool complete_clusters = true;
  std::vector<arma::mat> v_inverse_d_matrices;
  std::vector<arma::vec> sd_vectors;
  v_inverse_d_matrices.reserve(clusters_no);
  sd_vectors.reserve(clusters_no);
  for (arma::uword cluster_index = 0; cluster_index < clusters_no;
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      v_inverse_d_matrices.push_back(
        solve_chol_or_lu_mat(cluster.v_matrix, cluster.d_matrix));
      sd_vectors.push_back(arma::sqrt(arma::diagvec(cluster.v_matrix)));
      const arma::mat& v_inverse_d_matrix = v_inverse_d_matrices.back();
      const arma::vec u_vector = v_inverse_d_matrix.t() * cluster.residuals;
      information_matrix += v_inverse_d_matrix.t() * cluster.d_matrix;
      meat_matrix += u_vector * u_vector.t();
      if (cluster.residuals.n_elem != repeated_max) complete_clusters = false;
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_small_sample_covariances",
                                   cluster_index, cluster.id_value, e);
    }
  }
  symmetrize_if_close(information_matrix, "get_small_sample_covariances");
  const arma::mat naive_covariance = solve_chol_or_lu_mat(
    information_matrix, arma::eye(params_no, params_no),
    "get_small_sample_covariances");
  arma::mat naive_root;
  const arma::mat naive_symmetric =
    0.5 * (naive_covariance + naive_covariance.t());
  if (!arma::chol(naive_root, naive_symmetric)) {
    Rcpp::stop("get_small_sample_covariances: the naive covariance matrix is "
               "not positive definite");
  }
  std::vector<double> singular_ids;

  arma::mat meat_kc(params_no, params_no, arma::fill::zeros);
  arma::mat meat_md(params_no, params_no, arma::fill::zeros);
  arma::mat meat_fg(params_no, params_no, arma::fill::zeros);
  arma::mat meat_fz(params_no, params_no, arma::fill::zeros);
  arma::mat pooled_pan(repeated_max, repeated_max, arma::fill::zeros);
  arma::mat pooled_wl(repeated_max, repeated_max, arma::fill::zeros);
  arma::mat pooled_wb(repeated_max, repeated_max, arma::fill::zeros);
  for (arma::uword cluster_index = 0; cluster_index < clusters_no;
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::mat& d_matrix = cluster.d_matrix;
      const arma::mat& v_inverse_d_matrix = v_inverse_d_matrices[cluster_index];
      const arma::vec& residuals = cluster.residuals;
      const arma::uword cluster_size = residuals.n_elem;
      const arma::vec u_vector = v_inverse_d_matrix.t() * residuals;
      const bool leverage_one =
        get_max_cluster_leverage(naive_root, v_inverse_d_matrix, d_matrix) >
        1.0 - leverage_tolerance;
      arma::mat md_matrix;
      arma::vec md_residuals;
      arma::vec kc_residuals;
      if (leverage_one) {
        singular_ids.push_back(cluster.id_value);
      } else {
        const arma::mat h_matrix =
          d_matrix * naive_covariance * v_inverse_d_matrix.t();
        const arma::mat identity_matrix = arma::eye(cluster_size, cluster_size);
        const CholOrLuFactor md_factor(identity_matrix - h_matrix,
                                       "get_small_sample_covariances: I - H",
                                       false);
        md_matrix = md_factor.solve(identity_matrix);
        const arma::cx_mat kc_complex = arma::sqrtmat(md_matrix);
        if (arma::abs(arma::imag(kc_complex)).max() >
            1e-8 * std::max(1.0, arma::abs(arma::real(kc_complex)).max())) {
          throw std::runtime_error(
            "the matrix square root of (I - H)^{-1} is not real");
        }
        const arma::mat kc_matrix = arma::real(kc_complex);
        md_residuals = md_matrix * residuals;
        kc_residuals = kc_matrix * residuals;
        const arma::vec md_u_vector = v_inverse_d_matrix.t() * md_residuals;
        const arma::vec kc_u_vector = v_inverse_d_matrix.t() * kc_residuals;
        meat_md += md_u_vector * md_u_vector.t();
        meat_kc += kc_u_vector * kc_u_vector.t();
      }

      const arma::vec leverage_vector =
        arma::diagvec(v_inverse_d_matrix.t() * d_matrix * naive_covariance);
      const arma::vec fg_vector =
        1.0 / arma::sqrt(1.0 - arma::clamp(leverage_vector, -arma::datum::inf, 0.75));
      const arma::vec fg_u_vector = fg_vector % u_vector;
      meat_fg += fg_u_vector * fg_u_vector.t();

      if (!leverage_one) {
        const arma::mat inner_matrix =
          residuals * residuals.t() -
          d_matrix * naive_covariance * (meat_matrix - u_vector * u_vector.t()) *
            naive_covariance * d_matrix.t();
        const arma::mat fz_matrix = v_inverse_d_matrix.t() * md_matrix;
        meat_fz += fz_matrix * inner_matrix * fz_matrix.t();
      }

      if (complete_clusters) {
        const arma::vec pan_residuals = residuals / sd_vectors[cluster_index];
        pooled_pan += pan_residuals * pan_residuals.t();
        if (!leverage_one) {
          const arma::vec& sd_vector = sd_vectors[cluster_index];
          const arma::vec wl_residuals = md_residuals / sd_vector;
          const arma::vec wb_residuals = kc_residuals / sd_vector;
          pooled_wl += wl_residuals * wl_residuals.t();
          pooled_wb += wb_residuals * wb_residuals.t();
        }
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
    for (arma::uword cluster_index = 0; cluster_index < clusters_no;
         ++cluster_index) {
      arma::mat scaled_matrix = v_inverse_d_matrices[cluster_index];
      scaled_matrix.each_col() %= sd_vectors[cluster_index];
      meat_pan += scaled_matrix.t() * pooled_pan * scaled_matrix;
      meat_wl += scaled_matrix.t() * pooled_wl * scaled_matrix;
      meat_wb += scaled_matrix.t() * pooled_wb * scaled_matrix;
    }
  }

  if (!singular_ids.empty()) {
    meat_kc = na_matrix;
    meat_md = na_matrix;
    meat_fz = na_matrix;
    meat_wl = na_matrix;
    meat_wb = na_matrix;
  }
  auto sandwich = [&naive_covariance](const arma::mat& meat) {
    arma::mat covariance = naive_covariance * meat * naive_covariance;
    symmetrize_if_close(covariance, "get_small_sample_covariances");
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
    Rcpp::Named("westgate-burchett")   = sandwich(meat_wb),
    Rcpp::Named("singular_clusters")   = arma::vec(singular_ids)
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
  reset_lu_fallback_warning();
  check_cluster_inputs(id_vector,
                       {y_vector, repeated_vector, weights_vector, mu_vector,
                        eta_vector},
                       model_matrix);
  const auto clusters = clusters_from_sorted_id(id_vector);
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec delta_vector = link_derivative_1(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const arma::mat correlation_matrix =
    get_correlation_matrix_cc(correlation_structure, alpha_vector, repeated_max);
  std::vector<ClusterData> quantities;
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
      quantities.push_back({cluster.id_value,
                            d_matrix_i,
                            v_matrix_i,
                            s_vector.subvec(first_row, last_row)});
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
  reset_lu_fallback_warning();
  check_cluster_inputs(id_vector,
                       {y_vector, repeated_vector, weights_vector, mu_vector,
                        eta_vector},
                       model_matrix);
  const auto clusters = clusters_from_sorted_id(id_vector);
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec delta_vector = link_derivative_1(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  std::vector<ClusterData> quantities;
  quantities.reserve(clusters.size());
  for (arma::uword cluster_index = 0; cluster_index < clusters.size(); ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      arma::mat d_matrix_i = model_matrix.rows(first_row, last_row);
      d_matrix_i.each_col() %= delta_vector.subvec(first_row, last_row);
      const arma::vec odds_ratios_vector_i =
        get_subject_specific_odds_ratios_or(repeated_vector.subvec(first_row, last_row),
                                         repeated_max,
                                         alpha_vector);
      const arma::mat v_matrix_i =
        get_v_matrix_or(mu_vector.subvec(first_row, last_row),
                        odds_ratios_vector_i,
                        weights_vector.subvec(first_row, last_row));
      quantities.push_back({cluster.id_value,
                            d_matrix_i,
                            v_matrix_i,
                            s_vector.subvec(first_row, last_row)});
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_small_sample_covariances_or", cluster_index,
                                   id_vector[cluster.start], e);
    }
  }
  return get_small_sample_covariances(quantities, model_matrix.n_cols, repeated_max);
}
//==============================================================================
