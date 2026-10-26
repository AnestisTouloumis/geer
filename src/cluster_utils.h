#ifndef GEER_CLUSTER_UTILS_H
#define GEER_CLUSTER_UTILS_H

#include <RcppArmadillo.h>
#include "utils.h"
#include "family_codes.h"
#include <algorithm>
#include <vector>

struct Cluster {
  arma::uword first_row;
  arma::uword after_last_row;
  double id_value;
};

struct ClusterData {
  double id_value;
  arma::mat d_matrix;
  arma::mat v_matrix;
  arma::vec residuals;
};

struct ClusterWorkspace {
  ClusterData data;
  arma::vec odds_ratios;
};

arma::mat build_cluster_covariance_cc(
  const Cluster& cluster, const arma::vec& mu_vector,
  const arma::vec& repeated_vector, const arma::vec& weights_vector,
  FamilyCode family_code, double phi, const arma::mat& correlation_matrix);
ClusterWorkspace build_cluster_data_cc(
  const Cluster& cluster, const arma::mat& model_matrix,
  const arma::vec& mu_vector, const arma::vec& repeated_vector,
  const arma::vec& weights_vector, const arma::vec& delta_vector,
  const arma::vec& residuals, FamilyCode family_code, double phi,
  const arma::mat& correlation_matrix);
ClusterWorkspace build_cluster_data_or(
  const Cluster& cluster, const arma::mat& model_matrix,
  const arma::vec& mu_vector, const arma::vec& repeated_vector,
  const arma::vec& weights_vector, const arma::vec& delta_vector,
  const arma::vec& residuals, arma::uword repeated_max,
  const arma::vec& alpha_vector);
arma::mat build_cluster_covariance_or(
  const Cluster& cluster, const arma::vec& mu_vector,
  const arma::vec& repeated_vector, const arma::vec& weights_vector,
  arma::uword repeated_max, const arma::vec& alpha_vector);

inline std::vector<Cluster> clusters_from_sorted_id(const arma::vec& id) {
  const arma::uword n = id.n_elem;
  std::vector<Cluster> clusters;
  if (n == 0) return clusters;
  if (!id.is_finite()) {
    signal_input_error("'id' must contain only finite values");
  }
  arma::uword n_clusters = 1;
  for (arma::uword i = 1; i < n; ++i) {
    if (id[i] != id[i - 1]) ++n_clusters;
  }
  clusters.reserve(n_clusters);
  arma::uword first_row = 0;
  while (first_row < n) {
    arma::uword after_last_row = first_row + 1;
    while (after_last_row < n && id[after_last_row] == id[first_row]) {
      ++after_last_row;
    }
    clusters.push_back({first_row, after_last_row, id[first_row]});
    first_row = after_last_row;
  }
  if (clusters.size() > 1) {
    std::vector<double> cluster_ids;
    cluster_ids.reserve(clusters.size());
    for (const auto& cluster : clusters) {
      cluster_ids.push_back(id[cluster.first_row]);
    }
    std::sort(cluster_ids.begin(), cluster_ids.end());
    if (std::adjacent_find(cluster_ids.begin(), cluster_ids.end()) != cluster_ids.end()) {
      signal_input_error(
      "the observations of each cluster must be in adjacent"
      " rows; an 'id' value "
                 "occurs in non-adjacent rows");
    }
  }
  return clusters;
}

#endif
