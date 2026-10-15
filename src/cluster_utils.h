#ifndef GEER_CLUSTER_UTILS_H
#define GEER_CLUSTER_UTILS_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>

struct Cluster { arma::uword start, end; };  // [start, end)

inline std::vector<Cluster> clusters_from_sorted_id(const arma::vec& id) {
  const arma::uword n = id.n_elem;
  std::vector<Cluster> clusters;
  if (n == 0) return clusters;
  if (!id.is_finite()) {
    Rcpp::stop("'id' must contain only finite values");
  }
  arma::uword n_clusters = 1;
  for (arma::uword i = 1; i < n; ++i) {
    if (id[i] != id[i - 1]) ++n_clusters;
  }
  clusters.reserve(n_clusters);
  arma::uword s = 0;
  while (s < n) {
    arma::uword e = s + 1;
    while (e < n && id[e] == id[s]) ++e;
    clusters.push_back({s, e});
    s = e;
  }
  // The observations of a cluster must be contiguous: an id value that
  // reappears after another id would otherwise be silently split into several
  // clusters and give a wrong fit.
  if (clusters.size() > 1) {
    std::vector<double> cluster_ids;
    cluster_ids.reserve(clusters.size());
    for (const auto& cluster : clusters) cluster_ids.push_back(id[cluster.start]);
    std::sort(cluster_ids.begin(), cluster_ids.end());
    if (std::adjacent_find(cluster_ids.begin(), cluster_ids.end()) != cluster_ids.end()) {
      Rcpp::stop("the observations of each cluster must be in adjacent rows; an 'id' value "
                 "occurs in non-adjacent rows");
    }
  }
  return clusters;
}

#endif
