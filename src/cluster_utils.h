#ifndef CLUSTER_UTILS_H
#define CLUSTER_UTILS_H

#include <RcppArmadillo.h>
#include <vector>

struct Cluster { arma::uword start, end; };  // [start, end)

inline std::vector<Cluster> clusters_from_sorted_id(const arma::vec& id) {
  const arma::uword n = id.n_elem;
  std::vector<Cluster> clusters;
  if (n == 0) return clusters;
  if (!id.is_finite()) {
    Rcpp::stop("'id' must contain only finite values.");
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
  return clusters;
}

#endif
