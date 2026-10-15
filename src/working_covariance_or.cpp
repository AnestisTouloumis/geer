#include "working_covariance_or.h"
#include "utils.h"
#include <algorithm>
#include <cmath>


//============================ subject specific odds ratios ====================
arma::vec get_subject_specific_odds_ratios(const arma::vec& repeated_vector_i,
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

//============================ bivariate distribution ==========================
// Plackett distribution F(p, q; psi) = P(Y_1 = 1, Y_2 = 1) of two binary
// variables with marginal means p, q and odds ratio psi. With
//   f = 1 + (psi - 1) (p + q),  D = f^2 - 4 psi (psi - 1) p q,
// the textbook form F = (f - sqrt(D)) / (2 (psi - 1)) subtracts two nearly equal
// numbers (and is 0/0 at psi = 1). Because (f - sqrt(D)) (f + sqrt(D)) =
// 4 psi (psi - 1) p q, the algebraically identical F = 2 psi p q / (f + sqrt(D))
// is used whenever f > 0; it is smooth at psi = 1. For f <= 0 (only possible
// when psi < 1) f - sqrt(D) has no cancellation and the textbook form is used.
// Means within 1e-8 of 0 or 1 are treated as independent, as before.
bool is_plackett_boundary(const double row_prob, const double col_prob) {
  const double tol = 1e-8;
  return row_prob > 1.0 - tol || col_prob > 1.0 - tol ||
         row_prob < tol || col_prob < tol;
}

double get_bivariate_distribution(const double row_prob,
                                  const double col_prob,
                                  const double odds_ratio) {
  const double ans_independence = row_prob * col_prob;
  if (is_plackett_boundary(row_prob, col_prob)) {
    return ans_independence;
  }
  const double eps = odds_ratio - 1.0;
  const double f_value = 1.0 + eps * (row_prob + col_prob);
  const double root_value =
    std::sqrt(std::max(0.0, f_value * f_value - 4.0 * odds_ratio * eps * ans_independence));
  if (f_value > 0.0) {
    return 2.0 * odds_ratio * ans_independence / (f_value + root_value);
  }
  return (f_value - root_value) / (2.0 * eps);
}
//==============================================================================


//============================ bivariate covariance ============================
// Cov(Y_1, Y_2) = F - p q, evaluated without the cancellation in F - p q that
// makes the difference lose all relative accuracy as psi -> 1. With
// a = 2 psi - f one has F - p q = p q (a - sqrt(D)) / (f + sqrt(D)), and for
// a > 0, a - sqrt(D) = (a^2 - D) / (a + sqrt(D)) = 4 psi (psi - 1) (1 - p) (1 - q)
// / (a + sqrt(D)).
double get_bivariate_covariance(const double row_prob,
                                const double col_prob,
                                const double odds_ratio) {
  if (is_plackett_boundary(row_prob, col_prob)) {
    return 0.0;
  }
  const double ans_independence = row_prob * col_prob;
  const double eps = odds_ratio - 1.0;
  const double f_value = 1.0 + eps * (row_prob + col_prob);
  const double root_value =
    std::sqrt(std::max(0.0, f_value * f_value - 4.0 * odds_ratio * eps * ans_independence));
  if (f_value > 0.0) {
    const double a_value = 2.0 * odds_ratio - f_value;
    const double numerator = a_value > 0.0 ?
      4.0 * odds_ratio * eps * (1.0 - row_prob) * (1.0 - col_prob) /
        (a_value + root_value) :
      a_value - root_value;
    return ans_independence * numerator / (f_value + root_value);
  }
  return (f_value - root_value) / (2.0 * eps) - ans_independence;
}
//==============================================================================


}  // namespace

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
          get_bivariate_covariance(mu_vector[i], mu_vector[j], odds_ratios_vector[k]);
        ++k;
      }
    }
    result = arma::symmatu(result);
  }
  const arma::vec weights_sq_inv_vector = 1.0 / arma::sqrt(weights_vector);
  result %= (weights_sq_inv_vector * weights_sq_inv_vector.t());
  return result;
}
//==============================================================================


namespace {

//============================ first derivative wrt row probability ============
double get_bivariate_distribution_d_mu_row(const double row_prob,
                                        const double col_prob,
                                        const double odds_ratio) {
  const double f_value = 1.0 - (1.0 - odds_ratio) * (row_prob + col_prob);
  const double biv_dis =
    get_bivariate_distribution(row_prob, col_prob, odds_ratio);
  const double num =
    row_prob + col_prob - odds_ratio * (row_prob - col_prob) - 1.0;
  const double den =
    f_value - 2.0 * (odds_ratio - 1.0) * biv_dis;
  return 0.5 * (1.0 + num / den);
}
//==============================================================================


//============================ second derivative wrt row probability ===========
double get_bivariate_distribution_d2_mu_row(const double row_prob,
                                         const double col_prob,
                                         const double odds_ratio) {
  const double f_value = 1.0 - (1.0 - odds_ratio) * (row_prob + col_prob);
  const double biv_dis =
    get_bivariate_distribution(row_prob, col_prob, odds_ratio);
  const double num =
    2.0 * odds_ratio * (odds_ratio - 1.0) * col_prob * (col_prob - 1.0);
  const double den =
    std::pow(f_value - 2.0 * (odds_ratio - 1.0) * biv_dis, 3.0);
  return num / den;
}
//==============================================================================


//============================ second derivative wrt row-col probabilities =====
double get_bivariate_distribution_d2_mu_row_col(const double row_prob,
                                           const double col_prob,
                                           const double odds_ratio) {
  const double f_value = 1.0 - (1.0 - odds_ratio) * (row_prob + col_prob);
  const double biv_dis =
    get_bivariate_distribution(row_prob, col_prob, odds_ratio);
  const double num =
    (f_value - 2.0 * (odds_ratio - 1.0) * row_prob * col_prob) * odds_ratio;
  const double den =
    std::pow(f_value - 2.0 * (odds_ratio - 1.0) * biv_dis, 3.0);
  return num / den;
}
//==============================================================================


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
            upper_triangular_pair_index(smaller_index, larger_index, cluster_size);
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
  // No special case for cluster_size == 1: the -mu terms and the unit diagonal
  // below are then the only contributions, giving d V / d mu = 1 - 2 mu for
  // V = mu * (1 - mu), and the pair loop is simply not entered.
  for (arma::uword r = 0; r < cluster_size; ++r) {
    for (arma::uword i = 0; i < cluster_size; ++i) {
      result(r * cluster_size + i, i) = -mu_vector[r];
      result(r * cluster_size + i, r) -= mu_vector[i];

      if (i != r) {
        const arma::uword smaller_index = std::min(i, r);
        const arma::uword larger_index = std::max(i, r);
        const arma::uword k =
          upper_triangular_pair_index(smaller_index, larger_index, cluster_size);

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
  // kron(W^(-1/2), W^(-1/2)) is diagonal, with entry w_r^(-1/2) * w_i^(-1/2) at
  // row r * m + i, so the product is a row scaling.
  const arma::vec weights_sq_inv_vector = 1.0 / arma::sqrt(weights_vector);
  result.each_col() %=
    arma::vectorise(weights_sq_inv_vector * weights_sq_inv_vector.t());
  return result;
}
//==============================================================================

