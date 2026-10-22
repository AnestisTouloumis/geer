#include "nuisance_estimation_or.h"
#include "cluster_utils.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>


//============================ estimate marginalized odds ratio structure ======
// [[Rcpp::export]]
arma::vec get_marginalized_odds_ratios_or(const arma::vec& response_vector,
                                                 const arma::vec& id_vector,
                                                 const arma::vec& repeated_vector,
                                                 const arma::vec& weights_vector,
                                                 const double adding_constant,
                                                 const char* or_structure) {
  reset_lu_fallback_warning();
  const auto clusters = clusters_from_sorted_id(id_vector);
  const arma::uword sample_size = clusters.size();
  const arma::uword cluster_size_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::uword pairs_no = upper_triangular_pairs(cluster_size_max);
  arma::mat wide_responses_matrix(sample_size, cluster_size_max, arma::fill::zeros);
  arma::mat wide_weights_matrix(sample_size, cluster_size_max, arma::fill::zeros);
  for (arma::uword cluster_index = 0; cluster_index < sample_size; ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    const arma::uword first_row = cluster.start;
    const arma::uword last_row = cluster.end - 1;
    const arma::vec response_vector_i = response_vector.subvec(first_row, last_row);
    const arma::vec repeated_vector_i = repeated_vector.subvec(first_row, last_row);
    const arma::vec weights_vector_i = weights_vector.subvec(first_row, last_row);
    const arma::uword cluster_size_i = cluster.end - cluster.start;
    for (arma::uword j = 0; j < cluster_size_i; ++j) {
      const arma::uword time_index =
        static_cast<arma::uword>(repeated_vector_i[j]) - 1;
      wide_responses_matrix(cluster_index, time_index) = response_vector_i[j] + 1.0;
      wide_weights_matrix(cluster_index, time_index) = weights_vector_i[j];
    }
  }
  arma::vec counts(4 * pairs_no, arma::fill::zeros);
  for (arma::uword cluster_index = 0; cluster_index < sample_size;
       ++cluster_index) {
    for (arma::uword i = 0; i + 1 < cluster_size_max; ++i) {
      const double response_i = wide_responses_matrix(cluster_index, i);
      const double weight_i = wide_weights_matrix(cluster_index, i);
      for (arma::uword j = i + 1; j < cluster_size_max; ++j) {
        const double response_j = wide_responses_matrix(cluster_index, j);
        const double weight_j = wide_weights_matrix(cluster_index, j);
        if ((response_i == 1.0 || response_i == 2.0) &&
            (response_j == 1.0 || response_j == 2.0)) {
          const arma::uword pair_index =
            upper_triangular_pair_index(i, j, cluster_size_max);
          const arma::uword cell_index =
            static_cast<arma::uword>((response_i - 1.0) * 2.0 + (response_j - 1.0));
          counts[cell_index * pairs_no + pair_index] += std::min(weight_i, weight_j);
        }
      }
    }
  }
  std::vector<bool> informative(pairs_no, false);
  for (arma::uword i = 0; i < pairs_no; ++i) {
    informative[i] = (counts[i] + counts[i + pairs_no] +
                      counts[i + 2 * pairs_no] + counts[i + 3 * pairs_no]) > 0.0;
  }
  counts += adding_constant;
  arma::vec result(pairs_no);
  for (arma::uword i = 0; i < pairs_no; ++i) {
    result[i] =
      (counts[i] * counts[i + 3 * pairs_no]) /
      (counts[i + pairs_no] * counts[i + 2 * pairs_no]);
  }
  if (std::strcmp(or_structure, "exchangeable") == 0) {
    double log_sum = 0.0;
    arma::uword informative_no = 0;
    for (arma::uword i = 0; i < pairs_no; ++i) {
      if (informative[i]) {
        log_sum += std::log(result[i]);
        ++informative_no;
      }
    }
    if (informative_no > 0) {
      result.fill(std::exp(log_sum / static_cast<double>(informative_no)));
    }
  }
  return result;
}
//==============================================================================
