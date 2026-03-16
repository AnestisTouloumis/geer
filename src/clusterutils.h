#ifndef CLUSTERUTILS_H
#define CLUSTERUTILS_H

#include <RcppArmadillo.h>
#include <vector>

struct Cluster { arma::uword start, end; };  // [start, end)

inline std::vector<Cluster> clusters_from_sorted_id(const arma::vec& id) {
  const arma::uword n = id.n_elem;
  std::vector<Cluster> clusters;
  clusters.reserve(static_cast<std::size_t>(arma::max(id)));

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
