#include "nuisance_estimation_or.h"
#include "cluster_utils.h"
#include "utils.h"
#include <algorithm>


//============================ estimate marginalized odds ratio structure ======
// [[Rcpp::export]]
Rcpp::NumericVector get_marginalized_odds_ratios(const arma::vec& response_vector,
                                                 const arma::vec& id_vector,
                                                 const arma::vec& repeated_vector,
                                                 const arma::vec& weights_vector,
                                                 const double adding_constant,
                                                 const Rcpp::String& or_structure) {
  const auto clusters = clusters_from_sorted_id(id_vector);
  const arma::uword sample_size = clusters.size();
  const arma::uword cluster_size_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::uword pairs_no = upper_triangular_pairs(cluster_size_max);
  arma::mat wide_responses_matrix(sample_size, cluster_size_max, arma::fill::zeros);
  arma::mat wide_weights_matrix(sample_size, cluster_size_max, arma::fill::zeros);
  for (arma::uword cluster_idx = 0; cluster_idx < sample_size; ++cluster_idx) {
    const auto& cluster = clusters[cluster_idx];
    const arma::uword first_row = cluster.start;
    const arma::uword last_row = cluster.end - 1;
    const arma::vec response_vector_i = response_vector.subvec(first_row, last_row);
    const arma::vec repeated_vector_i = repeated_vector.subvec(first_row, last_row);
    const arma::vec weights_vector_i = weights_vector.subvec(first_row, last_row);
    const arma::uword cluster_size_i = cluster.end - cluster.start;
    for (arma::uword j = 0; j < cluster_size_i; ++j) {
      const arma::uword time_index =
        static_cast<arma::uword>(repeated_vector_i[j]) - 1;
      wide_responses_matrix(cluster_idx, time_index) = response_vector_i[j] + 1.0;
      wide_weights_matrix(cluster_idx, time_index) = weights_vector_i[j];
    }
  }
  arma::vec counts(4 * pairs_no, arma::fill::zeros);
  for (arma::uword l = 0; l < sample_size; ++l) {
    for (arma::uword i = 0; i + 1 < cluster_size_max; ++i) {
      const double response_i = wide_responses_matrix(l, i);
      const double weight_i = wide_weights_matrix(l, i);
      for (arma::uword j = i + 1; j < cluster_size_max; ++j) {
        const double response_j = wide_responses_matrix(l, j);
        const double weight_j = wide_weights_matrix(l, j);
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
  counts += adding_constant;
  Rcpp::NumericVector result(pairs_no);
  for (arma::uword i = 0; i < pairs_no; ++i) {
    result[i] =
      (counts[i] * counts[i + 3 * pairs_no]) /
        (counts[i + pairs_no] * counts[i + 2 * pairs_no]);
  }
  if (or_structure == "exchangeable") {
    result = Rcpp::exp(Rcpp::rep(Rcpp::mean(Rcpp::log(result)), pairs_no));
  }
  return result;
}
//==============================================================================
