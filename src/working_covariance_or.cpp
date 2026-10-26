#include "working_covariance_or.h"
#include "utils.h"
#include <algorithm>
#include <cmath>


//============================ subject specific odds ratios ====================
arma::vec get_cluster_odds_ratios_or(const arma::vec& repeated_vector_i,
                                           const arma::uword cluster_size_max,
                                           const arma::vec& odds_ratios_vector) {
  const arma::uword cluster_size_i = repeated_vector_i.n_elem;
  const arma::uword pairs_no_i = upper_triangular_pairs(cluster_size_i);
  arma::vec result(pairs_no_i, arma::fill::zeros);
  arma::uword k = 0;
  for (arma::uword i = 0; i + 1 < cluster_size_i; ++i) {
    const arma::uword row_index =
      static_cast<arma::uword>(repeated_vector_i[i]) - 1;
    for (arma::uword j = i + 1; j < cluster_size_i; ++j) {
      const arma::uword col_index =
        static_cast<arma::uword>(repeated_vector_i[j]) - 1;
      const arma::uword pos =
        upper_triangular_pair_index(row_index, col_index, cluster_size_max);

      result[k] = odds_ratios_vector[pos];
      ++k;
    }
  }
  return result;
}
//==============================================================================


namespace {

struct PlackettPair {
  double row_prob;
  double col_prob;
  double row_complement;
  double col_complement;
  double reciprocal;
  double departure;
  double root;
  double difference;
  bool complemented;
};

double compute_sum_minus_one(const double first, const double second) {
  const double larger = std::max(first, second);
  const double smaller = std::min(first, second);
  const double sum = larger + smaller;
  const double error = smaller - (sum - larger);
  return (sum - 1.0) + error;
}

PlackettPair compute_plackett_pair(const double row_prob,
                                  const double col_prob,
                                  const double odds_ratio) {
  if (!(R_FINITE(row_prob) && row_prob >= 0.0 && row_prob <= 1.0)) {
    signal_input_error("'row_prob' must lie between 0 and 1, inclusive");
  }
  if (!(R_FINITE(col_prob) && col_prob >= 0.0 && col_prob <= 1.0)) {
    signal_input_error("'col_prob' must lie between 0 and 1, inclusive");
  }
  if (!(R_FINITE(odds_ratio) && odds_ratio > 0.0)) {
    signal_input_error("'odds_ratio' must be a positive finite number");
  }
  const bool complemented = odds_ratio < 1.0;
  const double reciprocal = complemented ? odds_ratio : 1.0 / odds_ratio;
  const double departure = complemented ? 1.0 - odds_ratio :
    (odds_ratio - 1.0) / odds_ratio;
  const double col = complemented ? 1.0 - col_prob : col_prob;
  const double col_complement = complemented ? col_prob : 1.0 - col_prob;
  const double row_complement = 1.0 - row_prob;
  const double difference = complemented ?
    compute_sum_minus_one(row_prob, col_prob) : row_prob - col_prob;
  const double root = std::hypot(
    departure * difference,
    std::sqrt(reciprocal) * std::sqrt(
      2.0 * departure *
        (row_prob * col_complement + col * row_complement) + reciprocal));
  if (!(R_FINITE(root) && root > 0.0)) {
    Rcpp::stop("compute_plackett_pair failed: the discriminant is not "
               "a representable positive value");
  }
  return {row_prob, col, row_complement, col_complement,
          reciprocal, departure, root, difference, complemented};
}

double get_bivariate_covariance(const double row_prob,
                                const double col_prob,
                                const double odds_ratio) {
  const auto pair = compute_plackett_pair(row_prob, col_prob, odds_ratio);
  if (row_prob == 0.0 || row_prob == 1.0 ||
      col_prob == 0.0 || col_prob == 1.0) {
    return 0.0;
  }
  const double first = pair.reciprocal + pair.departure *
    (pair.row_prob + pair.col_prob) + pair.root;
  const double second = pair.reciprocal + pair.departure *
    (pair.row_complement + pair.col_complement) + pair.root;
  const double covariance = (pair.row_prob / first) * pair.col_prob *
    (4.0 * pair.departure * (pair.row_complement / second) *
     pair.col_complement);
  if (!R_FINITE(covariance)) {
    Rcpp::stop("get_bivariate_covariance failed: the covariance is not finite");
  }
  return pair.complemented ? -covariance : covariance;
}


}  // namespace

//============================ covariance for R callers ========================
// [[Rcpp::export]]
arma::vec get_plackett_covariance_or(const double row_prob,
                                    const double col_prob,
                                    const double odds_ratio) {
  arma::vec result(1);
  result[0] = get_bivariate_covariance(row_prob, col_prob, odds_ratio);
  return result;
}
//==============================================================================


//============================ v matrix ========================================
arma::mat get_v_matrix_or(const arma::vec& mu_vector,
                          const arma::vec& odds_ratios_vector,
                          const arma::vec& weights_vector) {
  const arma::uword cluster_size = mu_vector.n_elem;
  arma::mat result = arma::diagmat(mu_vector % (1.0 - mu_vector));
  if (cluster_size > 1) {
    arma::uword k = 0;
    for (arma::uword i = 0; i + 1 < cluster_size; ++i) {
      for (arma::uword j = i + 1; j < cluster_size; ++j) {
        result(i, j) =
          get_bivariate_covariance(mu_vector[i], mu_vector[j],
                                   odds_ratios_vector[k]);
        ++k;
      }
    }
    result = arma::symmatu(result);
  }
  const arma::vec weights_sq_inv_vector = 1.0 / arma::sqrt(weights_vector);
  result %= (weights_sq_inv_vector * weights_sq_inv_vector.t());
  check_working_covariance(result, "get_v_matrix_or");
  return result;
}
//==============================================================================


namespace {

double get_bivariate_distribution_d_mu_row(const double row_prob,
                                          const double col_prob,
                                          const double odds_ratio) {
  const auto pair = compute_plackett_pair(row_prob, col_prob, odds_ratio);
  const double numerator = -pair.departure * pair.difference +
    pair.reciprocal * (2.0 * pair.col_prob - 1.0);
  const double direction = pair.complemented ? -1.0 : 1.0;
  const double signed_numerator = direction * numerator;
  if (signed_numerator >= 0.0) {
    return 0.5 * (1.0 + signed_numerator / pair.root);
  }
  return 2.0 * (pair.reciprocal / pair.root) *
    (pair.col_prob / (pair.root - signed_numerator)) * pair.col_complement;
}

double get_bivariate_distribution_d2_mu_row(const double row_prob,
                                           const double col_prob,
                                           const double odds_ratio) {
  const auto pair = compute_plackett_pair(row_prob, col_prob, odds_ratio);
  const double derivative = -2.0 * pair.departure *
    (pair.reciprocal / pair.root) * (pair.col_prob / pair.root) *
    (pair.col_complement / pair.root);
  if (!R_FINITE(derivative)) {
    Rcpp::stop("get_bivariate_distribution_d2_mu_row failed: the derivative "
               "is not finite");
  }
  return pair.complemented ? -derivative : derivative;
}

double get_bivariate_distribution_d2_mu_row_col(const double row_prob,
                                               const double col_prob,
                                               const double odds_ratio) {
  const auto pair = compute_plackett_pair(row_prob, col_prob, odds_ratio);
  const double numerator = pair.reciprocal + pair.departure *
    (pair.row_prob * pair.col_complement +
     pair.col_prob * pair.row_complement);
  const double derivative = (pair.reciprocal / pair.root) *
    (numerator / pair.root) / pair.root;
  if (!R_FINITE(derivative)) {
    Rcpp::stop("get_bivariate_distribution_d2_mu_row_col failed: the "
               "derivative is not finite");
  }
  return derivative;
}


}  // namespace

//============================ derivatives g_matrix ============================
arma::mat get_g_matrix_or(const arma::vec& mu_vector,
                       const arma::vec& odds_ratios_vector) {
  const arma::uword cluster_size = mu_vector.n_elem;
  arma::mat result(cluster_size, cluster_size, arma::fill::zeros);

  if (cluster_size > 1) {
    for (arma::uword i = 0; i + 1 < cluster_size; ++i) {
      for (arma::uword j = i + 1; j < cluster_size; ++j) {
        const arma::uword k =
          upper_triangular_pair_index(i, j, cluster_size);

        result(i, j) =
          get_bivariate_distribution_d_mu_row(mu_vector[i],
                                           mu_vector[j],
                                           odds_ratios_vector[k]);
        result(j, i) =
          get_bivariate_distribution_d_mu_row(mu_vector[j],
                                           mu_vector[i],
                                           odds_ratios_vector[k]);
      }
    }
  }
  return result;
}
//==============================================================================


//============================ derivative of g matrix ==========================
arma::mat get_g_matrix_mu_or(const arma::vec& mu_vector,
                          const arma::vec& odds_ratios_vector) {
  const arma::uword cluster_size = mu_vector.n_elem;
  arma::mat result(cluster_size * cluster_size, cluster_size, arma::fill::zeros);

  if (cluster_size > 1) {
    for (arma::uword r = 0; r < cluster_size; ++r) {
      for (arma::uword i = 0; i < cluster_size; ++i) {
        if (i != r) {
          const arma::uword smaller_index = std::min(i, r);
          const arma::uword larger_index = std::max(i, r);
          const arma::uword k =
            upper_triangular_pair_index(smaller_index, larger_index,
                                      cluster_size);
          result(r * cluster_size + i, r) =
            get_bivariate_distribution_d2_mu_row_col(mu_vector[i],
                                                mu_vector[r],
                                                odds_ratios_vector[k]);
          result(r * cluster_size + i, i) =
            get_bivariate_distribution_d2_mu_row(mu_vector[i],
                                              mu_vector[r],
                                              odds_ratios_vector[k]);
        }
      }
    }
  }
  return result;
}
//==============================================================================


//============================ derivative of V_i wrt the means, dV/dmu =========
arma::mat get_v_matrix_mu_or(const arma::vec& mu_vector,
                             const arma::vec& odds_ratios_vector,
                             const arma::vec& weights_vector) {
  const arma::uword cluster_size = mu_vector.n_elem;
  arma::mat result(cluster_size * cluster_size, cluster_size, arma::fill::zeros);
  for (arma::uword r = 0; r < cluster_size; ++r) {
    for (arma::uword i = 0; i < cluster_size; ++i) {
      result(r * cluster_size + i, i) = -mu_vector[r];
      result(r * cluster_size + i, r) -= mu_vector[i];

      if (i != r) {
        const arma::uword smaller_index = std::min(i, r);
        const arma::uword larger_index = std::max(i, r);
        const arma::uword k =
          upper_triangular_pair_index(smaller_index, larger_index,
                                      cluster_size);

        result(r * cluster_size + i, i) +=
          get_bivariate_distribution_d_mu_row(mu_vector[i],
                                           mu_vector[r],
                                           odds_ratios_vector[k]);
        result(r * cluster_size + i, r) +=
          get_bivariate_distribution_d_mu_row(mu_vector[r],
                                           mu_vector[i],
                                           odds_ratios_vector[k]);
      }
    }
  }
  for (arma::uword j = 0; j < cluster_size; ++j) {
    result(j * cluster_size + j, j) += 1.0;
  }
  const arma::vec weights_sq_inv_vector = 1.0 / arma::sqrt(weights_vector);
  result.each_col() %=
    arma::vectorise(weights_sq_inv_vector * weights_sq_inv_vector.t());
  return result;
}
//==============================================================================

