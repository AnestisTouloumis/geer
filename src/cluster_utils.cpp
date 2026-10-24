#include "cluster_utils.h"
#include "working_covariance_cc.h"
#include "working_covariance_or.h"

namespace {
ClusterData build_cluster_data(const Cluster& cluster,
                                const arma::mat& model_matrix,
                                const arma::vec& delta_vector,
                                const arma::vec& residuals) {
  const arma::uword first_row = cluster.first_row;
  const arma::uword last_row = cluster.after_last_row - 1;
  ClusterData result;
  result.id_value = cluster.id_value;
  result.d_matrix = model_matrix.rows(first_row, last_row);
  result.d_matrix.each_col() %= delta_vector.subvec(first_row, last_row);
  result.residuals = residuals.subvec(first_row, last_row);
  return result;
}
}  // namespace

arma::mat build_cluster_covariance_cc(
  const Cluster& cluster, const arma::vec& mu_vector,
  const arma::vec& repeated_vector, const arma::vec& weights_vector,
  FamilyCode family_code, double phi, const arma::mat& correlation_matrix) {
  const arma::uword first_row = cluster.first_row;
  const arma::uword last_row = cluster.after_last_row - 1;
  return get_v_matrix_cc(
    family_code, mu_vector.subvec(first_row, last_row),
    repeated_vector.subvec(first_row, last_row), phi, correlation_matrix,
    weights_vector.subvec(first_row, last_row));
}

ClusterWorkspace build_cluster_data_cc(
  const Cluster& cluster, const arma::mat& model_matrix,
  const arma::vec& mu_vector, const arma::vec& repeated_vector,
  const arma::vec& weights_vector, const arma::vec& delta_vector,
  const arma::vec& residuals, FamilyCode family_code, double phi,
  const arma::mat& correlation_matrix) {
  ClusterWorkspace result;
  result.data = build_cluster_data(
    cluster, model_matrix, delta_vector, residuals);
  result.data.v_matrix = build_cluster_covariance_cc(
    cluster, mu_vector, repeated_vector, weights_vector,
    family_code, phi, correlation_matrix);
  return result;
}

ClusterWorkspace build_cluster_data_or(
  const Cluster& cluster, const arma::mat& model_matrix,
  const arma::vec& mu_vector, const arma::vec& repeated_vector,
  const arma::vec& weights_vector, const arma::vec& delta_vector,
  const arma::vec& residuals, arma::uword repeated_max,
  const arma::vec& alpha_vector) {
  const arma::uword first_row = cluster.first_row;
  const arma::uword last_row = cluster.after_last_row - 1;
  ClusterWorkspace result;
  result.data = build_cluster_data(
    cluster, model_matrix, delta_vector, residuals);
  result.odds_ratios = get_cluster_odds_ratios_or(
    repeated_vector.subvec(first_row, last_row), repeated_max, alpha_vector);
  result.data.v_matrix = get_v_matrix_or(
    mu_vector.subvec(first_row, last_row), result.odds_ratios,
    weights_vector.subvec(first_row, last_row));
  return result;
}

arma::mat build_cluster_covariance_or(
  const Cluster& cluster, const arma::vec& mu_vector,
  const arma::vec& repeated_vector, const arma::vec& weights_vector,
  arma::uword repeated_max, const arma::vec& alpha_vector) {
  const arma::uword first_row = cluster.first_row;
  const arma::uword last_row = cluster.after_last_row - 1;
  const arma::vec odds_ratios = get_cluster_odds_ratios_or(
    repeated_vector.subvec(first_row, last_row), repeated_max, alpha_vector);
  return get_v_matrix_or(mu_vector.subvec(first_row, last_row), odds_ratios,
                         weights_vector.subvec(first_row, last_row));
}
