#include "influence_measures.h"
#include "cluster_utils.h"
#include "link_functions.h"
#include "utils.h"
#include "variance_functions.h"
#include "working_covariance_cc.h"
#include "working_covariance_or.h"


namespace {
constexpr double leverage_tolerance = 1e-8;


//============================ per-cluster quantities ==========================
struct InfluenceCluster {
  double id_value;
  arma::mat d_matrix;
  arma::mat v_matrix;
  arma::vec residuals;
  arma::mat sensitivity_matrix;
};
//==============================================================================


//============================ influence quantities ============================
Rcpp::List get_influence_quantities(const std::vector<InfluenceCluster>& clusters,
                                    const arma::uword params_no,
                                    const arma::uword obs_no) {
  const arma::uword clusters_no = clusters.size();
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat sensitivity_matrix(params_no, params_no, arma::fill::zeros);
  std::vector<arma::mat> v_inverse_matrices;
  v_inverse_matrices.reserve(clusters_no);
  for (arma::uword cluster_index = 0; cluster_index < clusters_no; ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword cluster_size = cluster.residuals.n_elem;
      v_inverse_matrices.push_back(
        solve_chol_or_lu_mat(cluster.v_matrix,
                             arma::eye(cluster_size, cluster_size),
                             "get_influence_quantities"));
      information_matrix +=
        cluster.d_matrix.t() * v_inverse_matrices.back() * cluster.d_matrix;
      sensitivity_matrix += cluster.sensitivity_matrix;
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_influence_quantities", cluster_index,
                                   cluster.id_value, e);
    }
  }
  symmetrize_if_close(information_matrix, 1e-10);
  arma::mat naive_covariance = solve_chol_or_lu_mat(
    information_matrix, arma::eye(params_no, params_no),
    "get_influence_quantities");
  symmetrize_if_close(naive_covariance, 1e-10);
  arma::mat naive_root;
  if (!arma::chol(naive_root, naive_covariance)) {
    Rcpp::stop("get_influence_quantities: the naive covariance matrix is not "
               "positive definite");
  }

  arma::vec hat_vector(obs_no);
  arma::mat cluster_dfbeta(clusters_no, params_no);
  arma::mat observation_dfbeta(obs_no, params_no);
  arma::mat cluster_delta(clusters_no, params_no);
  arma::mat observation_delta(obs_no, params_no);
  arma::mat response_delta(obs_no, params_no);
  arma::uword first_row = 0;
  for (arma::uword cluster_index = 0; cluster_index < clusters_no; ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::mat& d_matrix = cluster.d_matrix;
      const arma::mat& v_inverse_matrix = v_inverse_matrices[cluster_index];
      const arma::vec& residuals = cluster.residuals;
      const arma::uword cluster_size = residuals.n_elem;
      const arma::uword last_row = first_row + cluster_size - 1;
      const arma::mat v_inverse_d_matrix = v_inverse_matrix * d_matrix;
      const arma::vec z_vector = v_inverse_matrix * residuals;
      const arma::mat scaled_matrix = v_inverse_d_matrix * naive_covariance;
      hat_vector.subvec(first_row, last_row) =
        arma::sum(d_matrix % scaled_matrix, 1);
      arma::mat cluster_leverage =
        naive_root * v_inverse_d_matrix.t() * d_matrix * naive_root.t();
      cluster_leverage = 0.5 * (cluster_leverage + cluster_leverage.t());
      if (arma::max(arma::eig_sym(cluster_leverage)) > 1.0 - leverage_tolerance) {
        cluster_dfbeta.row(cluster_index).fill(NA_REAL);
      } else {
        const arma::mat h_matrix = d_matrix * scaled_matrix.t();
        const arma::vec md_residuals = solve_chol_or_lu_vec(
          arma::eye(cluster_size, cluster_size) - h_matrix, residuals,
          "get_influence_quantities");
        cluster_dfbeta.row(cluster_index) =
          (scaled_matrix.t() * md_residuals).t();
      }
      const arma::vec scaled_leverage =
        arma::sum(v_inverse_d_matrix % scaled_matrix, 1);
      for (arma::uword j = 0; j < cluster_size; ++j) {
        const double precision = v_inverse_matrix(j, j);
        if (scaled_leverage[j] > (1.0 - leverage_tolerance) * precision) {
          observation_dfbeta.row(first_row + j).fill(NA_REAL);
        } else {
          observation_dfbeta.row(first_row + j) =
            z_vector[j] / (precision - scaled_leverage[j]) * scaled_matrix.row(j);
        }
      }
      cluster_delta.row(cluster_index) = (d_matrix.t() * z_vector).t();
      arma::mat observation_block = d_matrix;
      observation_block.each_col() %= z_vector;
      observation_delta.rows(first_row, last_row) = observation_block;
      arma::mat response_block = v_inverse_d_matrix;
      response_block.each_col() %= arma::sqrt(arma::diagvec(cluster.v_matrix));
      response_delta.rows(first_row, last_row) = response_block;
      first_row = last_row + 1;
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_influence_quantities", cluster_index,
                                   cluster.id_value, e);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("hat")                 = hat_vector,
    Rcpp::Named("cluster_dfbeta")      = cluster_dfbeta,
    Rcpp::Named("observation_dfbeta")  = observation_dfbeta,
    Rcpp::Named("sensitivity")         = sensitivity_matrix,
    Rcpp::Named("cluster_delta")       = cluster_delta,
    Rcpp::Named("observation_delta")   = observation_delta,
    Rcpp::Named("response_delta")      = response_delta
  );
}
//==============================================================================
}  // namespace


//============================ influence quantities -- cc ======================
// [[Rcpp::export]]
Rcpp::List get_influence_quantities_cc(const arma::vec& y_vector,
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
  const LinkCode link_code = parse_link(link);
  const FamilyCode family_code = parse_family(family);
  const auto clusters = clusters_from_sorted_id(id_vector);
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec delta_vector = link_derivative_1(link_code, eta_vector);
  const arma::vec delta_2_vector = link_derivative_2(link_code, eta_vector);
  const arma::vec variance_vector = variance_function(family_code, mu_vector);
  const arma::vec variance_1_vector = variance_derivative_1(family_code, mu_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
  std::vector<InfluenceCluster> quantities;
  quantities.reserve(clusters.size());
  for (arma::uword cluster_index = 0; cluster_index < clusters.size(); ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      const arma::mat x_matrix_i = model_matrix.rows(first_row, last_row);
      const arma::vec delta_i = delta_vector.subvec(first_row, last_row);
      const arma::vec residuals_i = s_vector.subvec(first_row, last_row);
      const arma::vec ratio_i = variance_1_vector.subvec(first_row, last_row) /
        variance_vector.subvec(first_row, last_row);
      arma::mat d_matrix_i = x_matrix_i;
      d_matrix_i.each_col() %= delta_i;
      const arma::mat v_matrix_i =
        get_v_matrix_cc(family_code,
                        mu_vector.subvec(first_row, last_row),
                        repeated_vector.subvec(first_row, last_row),
                        phi,
                        correlation_matrix,
                        weights_vector.subvec(first_row, last_row));
      const arma::vec z_vector_i =
        solve_chol_or_lu_vec(v_matrix_i, residuals_i, "get_influence_quantities_cc");
      const arma::vec a_vector = delta_2_vector.subvec(first_row, last_row) -
        0.5 * arma::square(delta_i) % ratio_i;
      const arma::vec b_vector = delta_i % (1.0 + 0.5 * residuals_i % ratio_i);
      arma::mat weighted_x_matrix = x_matrix_i;
      weighted_x_matrix.each_col() %= b_vector;
      arma::mat curvature_x_matrix = x_matrix_i;
      curvature_x_matrix.each_col() %= a_vector % z_vector_i;
      const arma::mat sensitivity_matrix_i =
        d_matrix_i.t() *
          solve_chol_or_lu_mat(v_matrix_i, weighted_x_matrix,
                               "get_influence_quantities_cc") -
        x_matrix_i.t() * curvature_x_matrix;
      quantities.push_back({id_vector[first_row],
                            d_matrix_i,
                            v_matrix_i,
                            residuals_i,
                            sensitivity_matrix_i});
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_influence_quantities_cc", cluster_index,
                                   id_vector[cluster.start], e);
    }
  }
  return get_influence_quantities(quantities, model_matrix.n_cols, y_vector.n_elem);
}
//==============================================================================


//============================ influence quantities -- or ======================
// [[Rcpp::export]]
Rcpp::List get_influence_quantities_or(const arma::vec& y_vector,
                                       const arma::mat& model_matrix,
                                       const arma::vec& id_vector,
                                       const arma::vec& repeated_vector,
                                       const arma::vec& weights_vector,
                                       const char* link,
                                       const arma::vec& mu_vector,
                                       const arma::vec& eta_vector,
                                       const arma::vec& alpha_vector) {
  const LinkCode link_code = parse_link(link);
  const auto clusters = clusters_from_sorted_id(id_vector);
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec delta_vector = link_derivative_1(link_code, eta_vector);
  const arma::vec delta_2_vector = link_derivative_2(link_code, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  std::vector<InfluenceCluster> quantities;
  quantities.reserve(clusters.size());
  for (arma::uword cluster_index = 0; cluster_index < clusters.size(); ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      const arma::uword cluster_size = last_row - first_row + 1;
      const arma::mat x_matrix_i = model_matrix.rows(first_row, last_row);
      const arma::vec mu_vector_i = mu_vector.subvec(first_row, last_row);
      const arma::vec weights_vector_i = weights_vector.subvec(first_row, last_row);
      const arma::vec residuals_i = s_vector.subvec(first_row, last_row);
      arma::mat d_matrix_i = x_matrix_i;
      d_matrix_i.each_col() %= delta_vector.subvec(first_row, last_row);
      const arma::vec odds_ratios_vector_i =
        get_subject_specific_odds_ratios(repeated_vector.subvec(first_row, last_row),
                                         repeated_max,
                                         alpha_vector);
      const arma::mat v_matrix_i =
        get_v_matrix_or(mu_vector_i, odds_ratios_vector_i, weights_vector_i);
      const arma::vec z_vector_i =
        solve_chol_or_lu_vec(v_matrix_i, residuals_i, "get_influence_quantities_or");
      const arma::mat v_matrix_mu_i =
        get_v_matrix_mu_or(mu_vector_i, odds_ratios_vector_i, weights_vector_i);
      arma::mat covariance_change_matrix(cluster_size, cluster_size);
      for (arma::uword m = 0; m < cluster_size; ++m) {
        covariance_change_matrix.col(m) =
          arma::reshape(v_matrix_mu_i.col(m), cluster_size, cluster_size) *
          z_vector_i;
      }
      arma::mat curvature_x_matrix = x_matrix_i;
      curvature_x_matrix.each_col() %=
        delta_2_vector.subvec(first_row, last_row) % z_vector_i;
      const arma::mat sensitivity_matrix_i =
        d_matrix_i.t() *
          solve_chol_or_lu_mat(v_matrix_i,
                               d_matrix_i + covariance_change_matrix * d_matrix_i,
                               "get_influence_quantities_or") -
        x_matrix_i.t() * curvature_x_matrix;
      quantities.push_back({id_vector[first_row],
                            d_matrix_i,
                            v_matrix_i,
                            residuals_i,
                            sensitivity_matrix_i});
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_influence_quantities_or", cluster_index,
                                   id_vector[cluster.start], e);
    }
  }
  return get_influence_quantities(quantities, model_matrix.n_cols, y_vector.n_elem);
}
//==============================================================================
