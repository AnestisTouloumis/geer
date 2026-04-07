#include "nuisance_quantities_or.h"
#include "clusterutils.h"
#include <algorithm>
#include <cmath>

namespace {
inline arma::uword upper_triangular_pairs(const arma::uword n) {
  return n * (n - 1) / 2;
}
  inline arma::uword upper_triangular_pair_index(const arma::uword i,
                                                 const arma::uword j,
                                                 const arma::uword n) {
    return i * (2 * n - i - 1) / 2 + (j - i - 1);
  }
}


//============================ estimate marginalized odds ratio structure ======
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector get_marginalized_odds_ratios(const arma::vec& response_vector,
                                                 const arma::vec& id_vector,
                                                 const arma::vec& repeated_vector,
                                                 const arma::vec& weights_vector,
                                                 const double& adding_constant,
                                                 const Rcpp::String& or_structure) {
  const auto clusters = clusters_from_sorted_id(id_vector);
  const arma::uword sample_size = clusters.size();
  const arma::uword cluster_size_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::uword pairs_no = upper_triangular_pairs(cluster_size_max);
  arma::mat wide_responses_matrix(sample_size, cluster_size_max, arma::fill::zeros);
  arma::mat wide_weights_matrix(sample_size, cluster_size_max, arma::fill::zeros);
  for (arma::uword cluster_idx = 0; cluster_idx < sample_size; ++cluster_idx) {
    const auto& cl = clusters[cluster_idx];
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::vec response_vector_i = response_vector.subvec(a, b);
    const arma::vec repeated_vector_i = repeated_vector.subvec(a, b);
    const arma::vec weights_vector_i = weights_vector.subvec(a, b);
    const arma::uword cluster_size_i = cl.end - cl.start;
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
  Rcpp::NumericVector ans(pairs_no);
  for (arma::uword i = 0; i < pairs_no; ++i) {
    ans[i] =
      (counts[i] * counts[i + 3 * pairs_no]) /
        (counts[i + pairs_no] * counts[i + 2 * pairs_no]);
  }
  if (or_structure == "exchangeable") {
    ans = Rcpp::exp(Rcpp::rep(Rcpp::mean(Rcpp::log(ans)), pairs_no));
  }
  return ans;
}
//==============================================================================


//============================ subject specific odds ratios ====================
// [[Rcpp::export]]
arma::vec get_subject_specific_odds_ratios(const arma::vec& repeated_vector_i,
                                           const int& cluster_size_max,
                                           const arma::vec& odds_ratios_vector) {
  const arma::uword cluster_size_i = repeated_vector_i.n_elem;
  const arma::uword pairs_no_i = upper_triangular_pairs(cluster_size_i);
  const arma::uword cluster_size_max_u =
    static_cast<arma::uword>(cluster_size_max);
  arma::vec ans(pairs_no_i, arma::fill::zeros);
  arma::uword k = 0;
  for (arma::uword i = 0; i + 1 < cluster_size_i; ++i) {
    const arma::uword row_index =
      static_cast<arma::uword>(repeated_vector_i[i]) - 1;
    for (arma::uword j = i + 1; j < cluster_size_i; ++j) {
      const arma::uword col_index =
        static_cast<arma::uword>(repeated_vector_i[j]) - 1;
      const arma::uword pos =
        upper_triangular_pair_index(row_index, col_index, cluster_size_max_u);

      ans[k] = odds_ratios_vector[pos];
      ++k;
    }
  }
  return ans;
}
//==============================================================================


//============================ bivariate distribution ==========================
// [[Rcpp::export]]
double get_bivariate_distribution(const double& row_prob,
                                  const double& col_prob,
                                  const double& odds_ratio) {
  const double ans_independence = row_prob * col_prob;
  const double tol = 1e-6;
  if (row_prob > 1.0 - tol || col_prob > 1.0 - tol) {
    return ans_independence;
  }
  if (row_prob < tol || col_prob < tol) {
    return ans_independence;
  }
  if (odds_ratio == 1.0) {
    return ans_independence;
  }
  const double f_value = 1.0 - (1.0 - odds_ratio) * (row_prob + col_prob);
  const double root_value = std::max(
    0.0,
    std::pow(f_value, 2.0) -
      4.0 * odds_ratio * (odds_ratio - 1.0) * ans_independence
  );
  return (f_value - std::sqrt(root_value)) / (2.0 * (odds_ratio - 1.0));
}
//==============================================================================


//============================ v matrix ========================================
// [[Rcpp::export]]
arma::mat get_v_matrix_or(const arma::vec& mu_vector,
                          const arma::vec& odds_ratios_vector,
                          const arma::vec& weights_vector) {
  const arma::uword cluster_size = mu_vector.n_elem;
  arma::mat ans = arma::diagmat(mu_vector % (1.0 - mu_vector));
  if (cluster_size > 1) {
    arma::uword k = 0;
    for (arma::uword i = 0; i + 1 < cluster_size; ++i) {
      for (arma::uword j = i + 1; j < cluster_size; ++j) {
        ans(i, j) =
          get_bivariate_distribution(mu_vector[i],
                                     mu_vector[j],
                                              odds_ratios_vector[k]) -
                                                mu_vector[i] * mu_vector[j];
        ++k;
      }
    }
    ans = arma::symmatu(ans);
  }
  const arma::vec weights_sq_inv_vector = 1.0 / arma::sqrt(weights_vector);
  ans %= (weights_sq_inv_vector * weights_sq_inv_vector.t());
  return ans;
}
//==============================================================================


//============================ inverse weight matrix ===========================
// [[Rcpp::export]]
arma::mat get_v_matrix_inverse_or(const arma::vec& mu_vector,
                                  const arma::vec& odds_ratios_vector,
                                  const arma::vec& weights_vector) {
  const arma::mat v_matrix =
    get_v_matrix_or(mu_vector, odds_ratios_vector, weights_vector);

  arma::mat ans;
  const bool ok = arma::inv(ans, v_matrix, arma::inv_opts::allow_approx);
  if (!ok) {
    Rcpp::stop("get_v_matrix_inverse_or: matrix inversion failed.");
  }
  return ans;
}
//==============================================================================


//============================ first derivative wrt row probability ============
// [[Rcpp::export]]
double get_bivariate_distribution_murow(const double& row_prob,
                                        const double& col_prob,
                                        const double& odds_ratio) {
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
// [[Rcpp::export]]
double get_bivariate_distribution_murow2(const double& row_prob,
                                         const double& col_prob,
                                         const double& odds_ratio) {
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
// [[Rcpp::export]]
double get_bivariate_distribution_murowcol(const double& row_prob,
                                           const double& col_prob,
                                           const double& odds_ratio) {
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


//============================ derivatives g_matrix ============================
// [[Rcpp::export]]
arma::mat get_g_matrix(const arma::vec& mu_vector,
                       const arma::vec& odds_ratios_vector) {
  const arma::uword cluster_size = mu_vector.n_elem;
  arma::mat ans(cluster_size, cluster_size, arma::fill::zeros);

  if (cluster_size > 1) {
    for (arma::uword i = 0; i + 1 < cluster_size; ++i) {
      for (arma::uword j = i + 1; j < cluster_size; ++j) {
        const arma::uword k =
          upper_triangular_pair_index(i, j, cluster_size);

        ans(i, j) =
          get_bivariate_distribution_murow(mu_vector[i],
                                           mu_vector[j],
                                                    odds_ratios_vector[k]);
        ans(j, i) =
          get_bivariate_distribution_murow(mu_vector[j],
                                           mu_vector[i],
                                                    odds_ratios_vector[k]);
      }
    }
  }
  return ans;
}
//==============================================================================


//============================ joint probability distribution derivative =======
// [[Rcpp::export]]
arma::mat get_bivariate_distribution_mu(const arma::vec& mu_vector,
                                        const arma::vec& odds_ratios_vector) {
  const arma::uword cluster_size = mu_vector.n_elem;
  arma::mat ans(cluster_size * cluster_size, cluster_size, arma::fill::zeros);
  if (cluster_size > 1) {
    for (arma::uword r = 0; r < cluster_size; ++r) {
      for (arma::uword i = 0; i < cluster_size; ++i) {
        if (i != r) {
          const arma::uword a = std::min(i, r);
          const arma::uword b = std::max(i, r);
          const arma::uword k =
            upper_triangular_pair_index(a, b, cluster_size);
          ans(r * cluster_size + i, i) =
            get_bivariate_distribution_murow(mu_vector[i],
                                             mu_vector[r],
                                                      odds_ratios_vector[k]);

          ans(r * cluster_size + i, r) +=
            get_bivariate_distribution_murow(mu_vector[r],
                                             mu_vector[i],
                                                      odds_ratios_vector[k]);
        }
      }
    }
  }
  return ans;
}
//==============================================================================


//============================ derivative of g matrix ==========================
// [[Rcpp::export]]
arma::mat get_g_matrix_mu(const arma::vec& mu_vector,
                          const arma::vec& odds_ratios_vector) {
  const arma::uword cluster_size = mu_vector.n_elem;
  arma::mat ans(cluster_size * cluster_size, cluster_size, arma::fill::zeros);

  if (cluster_size > 1) {
    for (arma::uword r = 0; r < cluster_size; ++r) {
      for (arma::uword i = 0; i < cluster_size; ++i) {
        if (i != r) {
          const arma::uword a = std::min(i, r);
          const arma::uword b = std::max(i, r);
          const arma::uword k =
            upper_triangular_pair_index(a, b, cluster_size);
          ans(r * cluster_size + i, r) =
            get_bivariate_distribution_murowcol(mu_vector[i],
                                                mu_vector[r],
                                                         odds_ratios_vector[k]);
          ans(r * cluster_size + i, i) =
            get_bivariate_distribution_murow2(mu_vector[i],
                                              mu_vector[r],
                                                       odds_ratios_vector[k]);
        }
      }
    }
  }
  return ans;
}
//==============================================================================


//============================ derivative of weight matrix =====================
// [[Rcpp::export]]
arma::mat get_v_matrix_mu_or(const arma::vec& mu_vector,
                             const arma::vec& odds_ratios_vector,
                             const arma::vec& weights_vector) {
  const arma::uword cluster_size = mu_vector.n_elem;
  arma::mat ans(cluster_size * cluster_size, cluster_size, arma::fill::zeros);
  if (cluster_size > 1) {
    for (arma::uword r = 0; r < cluster_size; ++r) {
      for (arma::uword i = 0; i < cluster_size; ++i) {
        ans(r * cluster_size + i, i) = -mu_vector[r];
        ans(r * cluster_size + i, r) -= mu_vector[i];

        if (i != r) {
          const arma::uword a = std::min(i, r);
          const arma::uword b = std::max(i, r);
          const arma::uword k =
            upper_triangular_pair_index(a, b, cluster_size);

          ans(r * cluster_size + i, i) +=
            get_bivariate_distribution_murow(mu_vector[i],
                                             mu_vector[r],
                                                      odds_ratios_vector[k]);
          ans(r * cluster_size + i, r) +=
            get_bivariate_distribution_murow(mu_vector[r],
                                             mu_vector[i],
                                                      odds_ratios_vector[k]);
        }
      }
    }
    for (arma::uword j = 0; j < cluster_size; ++j) {
      ans(j * cluster_size + j, j) += 1.0;
    }
  }
  const arma::mat weights_sq_inverse_matrix =
    arma::diagmat(1.0 / arma::sqrt(weights_vector));
  return arma::kron(weights_sq_inverse_matrix, weights_sq_inverse_matrix) * ans;
}
//==============================================================================
