#include "nuisance_quantities_cc.h"

#include "clusterutils.h"
#include "utils.h"
#include "variance_functions.h"

#include <cfloat>
#include <cmath>

namespace {

inline bool is_contiguous_1based(const arma::vec& r) {
  if (r.n_elem <= 1) {
    return true;
  }

  for (arma::uword i = 1; i < r.n_elem; ++i) {
    if (r[i] != r[i - 1] + 1.0) {
      return false;
    }
  }

  return true;
}

  inline void clamp_alpha_vector(arma::vec& alpha_vector) {
    const double upper_bound = 1.0 - 10.0 * DBL_EPSILON;
    const double lower_bound = 10.0 * DBL_EPSILON - 1.0;

    for (arma::uword i = 0; i < alpha_vector.n_elem; ++i) {
      if (alpha_vector[i] >= upper_bound) {
        alpha_vector[i] = upper_bound;
      }
      if (alpha_vector[i] <= lower_bound) {
        alpha_vector[i] = lower_bound;
      }
    }
  }

}  // namespace


//============================ pearson residuals ===============================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec get_pearson_residuals(const char* family,
                                const arma::vec& y_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& weights_vector) {
  const arma::vec variance_vector = variance(family, mu_vector);
  const arma::vec scale = arma::sqrt(weights_vector / variance_vector);
  return (y_vector - mu_vector) % scale;
}
//==============================================================================


//============================ phi hat =========================================
// [[Rcpp::export]]
double get_phi_hat(const arma::vec& pearson_residuals_vector,
                   const int& params_no) {
  const double n = static_cast<double>(pearson_residuals_vector.n_elem);
  const double denominator = n - static_cast<double>(params_no);

  if (denominator <= 0.0) {
    Rcpp::stop("get_phi_hat: non-positive denominator.");
  }

  double ans =
    arma::accu(arma::square(pearson_residuals_vector)) / denominator;

  if (ans < DBL_EPSILON) {
    ans = 10.0 * DBL_EPSILON;
  }

  return ans;
}
//==============================================================================


//============================ exchangeable alpha hat ==========================
// [[Rcpp::export]]
double alpha_hat_exchangeable(const arma::vec& pearson_residuals_vector,
                              const arma::vec& id_vector,
                              const double& phi,
                              const int& params_no) {
  double num = 0.0;
  double den = 0.0;

  const auto clusters = clusters_from_sorted_id(id_vector);

  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::vec pearson_residuals_vector_i =
      pearson_residuals_vector.subvec(a, b);
    const arma::uword cluster_size_i = cl.end - cl.start;

    for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
      for (arma::uword k = j + 1; k < cluster_size_i; ++k) {
        num += pearson_residuals_vector_i[j] * pearson_residuals_vector_i[k];
      }
    }

    den += static_cast<double>(cluster_size_i) *
      static_cast<double>(cluster_size_i - 1) * 0.5;
  }

  return num / ((den - static_cast<double>(params_no)) * phi);
}
//==============================================================================


//============================ ar1 alpha hat ===================================
// [[Rcpp::export]]
double alpha_hat_ar1(const arma::vec& pearson_residuals_vector,
                     const arma::vec& id_vector,
                     const arma::vec& repeated_vector,
                     const double& phi,
                     const int& params_no) {
  double num = 0.0;
  double den = 0.0;

  const auto clusters = clusters_from_sorted_id(id_vector);

  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::vec pearson_residuals_vector_i =
      pearson_residuals_vector.subvec(a, b);
    const arma::vec repeated_vector_i = repeated_vector.subvec(a, b);
    const arma::uword cluster_size_i = cl.end - cl.start;

    for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
      if (repeated_vector_i[j + 1] - repeated_vector_i[j] == 1.0) {
        num += pearson_residuals_vector_i[j] * pearson_residuals_vector_i[j + 1];
        den += 1.0;
      }
    }
  }

  return num / ((den - static_cast<double>(params_no)) * phi);
}
//==============================================================================


//============================ unstructured alpha hat ==========================
// [[Rcpp::export]]
arma::vec alpha_hat_unstructured(const arma::vec& pearson_residuals_vector,
                                 const arma::vec& id_vector,
                                 const arma::vec& repeated_vector,
                                 const double& phi,
                                 const int& params_no) {
  const arma::uword time_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::uword time_pairs = time_max * (time_max - 1) / 2;

  arma::vec num(time_pairs, arma::fill::zeros);
  arma::vec den(time_pairs, arma::fill::zeros);

  const auto clusters = clusters_from_sorted_id(id_vector);

  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::vec repeated_vector_i = repeated_vector.subvec(a, b);
    const arma::vec pearson_residuals_vector_i =
      pearson_residuals_vector.subvec(a, b);
    const arma::uword cluster_size_i = cl.end - cl.start;

    for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
      const arma::uword index_j =
        static_cast<arma::uword>(repeated_vector_i[j]);
      const arma::uword combn_j = index_j * (index_j + 1) / 2;

      for (arma::uword k = j + 1; k < cluster_size_i; ++k) {
        const arma::uword index_k =
          static_cast<arma::uword>(repeated_vector_i[k]);
        const arma::uword time_index =
          time_max * (index_j - 1) + index_k - combn_j;

        num[time_index - 1] +=
          pearson_residuals_vector_i[j] * pearson_residuals_vector_i[k];
        den[time_index - 1] += 1.0;
      }
    }
  }

  return num / ((den - static_cast<double>(params_no)) * phi);
}
//==============================================================================


//============================ m-dependent alpha hat ===========================
// [[Rcpp::export]]
arma::vec alpha_hat_mdependent(const arma::vec& pearson_residuals_vector,
                               const arma::vec& id_vector,
                               const arma::vec& repeated_vector,
                               const double& phi,
                               const int& params_no,
                               const int& mdependence) {
  arma::vec num(mdependence, arma::fill::zeros);
  arma::vec den(mdependence, arma::fill::zeros);

  const auto clusters = clusters_from_sorted_id(id_vector);

  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::vec repeated_vector_i = repeated_vector.subvec(a, b);
    const arma::vec pearson_residuals_vector_i =
      pearson_residuals_vector.subvec(a, b);
    const arma::uword cluster_size_i = cl.end - cl.start;

    for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
      const arma::uword index_j =
        static_cast<arma::uword>(repeated_vector_i[j]);

      for (arma::uword k = j + 1; k < cluster_size_i; ++k) {
        const arma::uword index_k =
          static_cast<arma::uword>(repeated_vector_i[k]);
        const arma::uword diff_int = index_k - index_j;

        if (diff_int < static_cast<arma::uword>(mdependence + 1)) {
          num[diff_int - 1] +=
            pearson_residuals_vector_i[j] * pearson_residuals_vector_i[k];
          den[diff_int - 1] += 1.0;
        }
      }
    }
  }

  return num / ((den - static_cast<double>(params_no)) * phi);
}
//==============================================================================


//============================ toeplitz alpha hat ==============================
// [[Rcpp::export]]
arma::vec alpha_hat_toeplitz(const arma::vec& pearson_residuals_vector,
                             const arma::vec& id_vector,
                             const arma::vec& repeated_vector,
                             const double& phi,
                             const int& params_no) {
  const arma::uword time_max =
    static_cast<arma::uword>(arma::max(repeated_vector));

  arma::vec num(time_max - 1, arma::fill::zeros);
  arma::vec den(time_max - 1, arma::fill::zeros);

  const auto clusters = clusters_from_sorted_id(id_vector);

  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::vec repeated_vector_i = repeated_vector.subvec(a, b);
    const arma::vec pearson_residuals_vector_i =
      pearson_residuals_vector.subvec(a, b);
    const arma::uword cluster_size_i = cl.end - cl.start;

    for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
      const arma::uword index_j =
        static_cast<arma::uword>(repeated_vector_i[j]);

      for (arma::uword k = j + 1; k < cluster_size_i; ++k) {
        const arma::uword index_k =
          static_cast<arma::uword>(repeated_vector_i[k]);
        const arma::uword lag = index_k - index_j;

        num[lag - 1] +=
          pearson_residuals_vector_i[j] * pearson_residuals_vector_i[k];
        den[lag - 1] += 1.0;
      }
    }
  }

  return num / ((den - static_cast<double>(params_no)) * phi);
}
//==============================================================================


//============================ alpha hat =======================================
// [[Rcpp::export]]
arma::vec get_alpha_hat(const char* correlation_structure,
                        const arma::vec& pearson_residuals_vector,
                        const arma::vec& id_vector,
                        const arma::vec& repeated_vector,
                        const double& phi,
                        const int& params_no,
                        const int& mdependence) {
  arma::vec ans;

  if (std::strcmp(correlation_structure, "independence") == 0) {
    ans.reset();
  } else if (std::strcmp(correlation_structure, "exchangeable") == 0) {
    ans.set_size(1);
    ans[0] = alpha_hat_exchangeable(pearson_residuals_vector,
                                    id_vector,
                                    phi,
                                    params_no);
  } else if (std::strcmp(correlation_structure, "ar1") == 0) {
    ans.set_size(1);
    ans[0] = alpha_hat_ar1(pearson_residuals_vector,
                           id_vector,
                           repeated_vector,
                           phi,
                           params_no);
  } else if (std::strcmp(correlation_structure, "m-dependent") == 0) {
    ans = alpha_hat_mdependent(pearson_residuals_vector,
                               id_vector,
                               repeated_vector,
                               phi,
                               params_no,
                               mdependence);
  } else if (std::strcmp(correlation_structure, "unstructured") == 0) {
    ans = alpha_hat_unstructured(pearson_residuals_vector,
                                 id_vector,
                                 repeated_vector,
                                 phi,
                                 params_no);
  } else if (std::strcmp(correlation_structure, "toeplitz") == 0) {
    ans = alpha_hat_toeplitz(pearson_residuals_vector,
                             id_vector,
                             repeated_vector,
                             phi,
                             params_no);
  } else {
    Rcpp::stop("Unsupported correlation structure.");
  }

  clamp_alpha_vector(ans);
  return ans;
}
//==============================================================================


//============================ independence ====================================
// [[Rcpp::export]]
arma::mat correlation_independence(const arma::uword dimension) {
  return arma::eye(dimension, dimension);
}
//==============================================================================


//============================ exchangeable ====================================
// [[Rcpp::export]]
arma::mat correlation_exchangeable(const arma::vec& alpha_vector,
                                   const arma::uword dimension) {
  arma::mat ans(dimension, dimension);
  ans.fill(alpha_vector[0]);
  ans.diag().fill(1.0);
  return ans;
}
//==============================================================================


//============================ ar1 =============================================
// [[Rcpp::export]]
arma::mat correlation_ar1(const arma::vec& alpha_vector,
                          const arma::uword dimension) {
  arma::vec ans(dimension, arma::fill::zeros);
  ans[0] = 1.0;

  for (arma::uword i = 1; i < dimension; ++i) {
    ans[i] = ans[i - 1] * alpha_vector[0];
  }

  return arma::toeplitz(ans);
}
//==============================================================================


//============================ m-dependent =====================================
// [[Rcpp::export]]
arma::mat correlation_mdependent(const arma::vec& alpha_vector,
                                 const arma::uword dimension) {
  const arma::uword k = alpha_vector.n_elem;

  if (k + 1 > dimension) {
    Rcpp::stop("correlation_mdependent: alpha_vector is too long.");
  }

  arma::vec toeplitz_vector(dimension, arma::fill::zeros);
  toeplitz_vector[0] = 1.0;
  if (k > 0) {
    toeplitz_vector.subvec(1, k) = alpha_vector;
  }

  return arma::toeplitz(toeplitz_vector);
}
//==============================================================================


//============================ toeplitz ========================================
// [[Rcpp::export]]
arma::mat correlation_toeplitz(const arma::vec& alpha_vector) {
  arma::vec toeplitz_vector(alpha_vector.n_elem + 1, arma::fill::zeros);
  toeplitz_vector[0] = 1.0;
  if (!alpha_vector.is_empty()) {
    toeplitz_vector.subvec(1, alpha_vector.n_elem) = alpha_vector;
  }

  return arma::toeplitz(toeplitz_vector);
}
//==============================================================================


//============================ unstructured ====================================
// [[Rcpp::export]]
arma::mat correlation_unstructured(const arma::vec& alpha_vector,
                                   const arma::uword dimension) {
  arma::mat ans_lt = arma::eye(dimension, dimension);
  const arma::uvec lt_indices = arma::trimatl_ind(arma::size(ans_lt), -1);
  ans_lt.elem(lt_indices) = alpha_vector;
  return arma::symmatl(ans_lt);
}
//==============================================================================


//============================ correlation matrix given rho vector =============
// [[Rcpp::export]]
arma::mat get_correlation_matrix(const char* correlation_structure,
                                 const arma::vec& alpha_vector,
                                 const arma::uword dimension) {
  if (std::strcmp(correlation_structure, "independence") == 0) {
    return correlation_independence(dimension);
  }
  if (std::strcmp(correlation_structure, "exchangeable") == 0) {
    return correlation_exchangeable(alpha_vector, dimension);
  }
  if (std::strcmp(correlation_structure, "ar1") == 0) {
    return correlation_ar1(alpha_vector, dimension);
  }
  if (std::strcmp(correlation_structure, "m-dependent") == 0) {
    return correlation_mdependent(alpha_vector, dimension);
  }
  if (std::strcmp(correlation_structure, "unstructured") == 0) {
    return correlation_unstructured(alpha_vector, dimension);
  }
  if (std::strcmp(correlation_structure, "toeplitz") == 0) {
    return correlation_toeplitz(alpha_vector);
  }
  if (std::strcmp(correlation_structure, "fixed") == 0) {
    return correlation_unstructured(alpha_vector, dimension);
  }

  Rcpp::stop("Unsupported correlation structure.");
}
//==============================================================================


//============================ subject-specific weight matrix ==================
// [[Rcpp::export]]
arma::mat get_v_matrix_cc(const char* family,
                          const arma::vec& mu_vector,
                          const arma::vec& repeated_vector,
                          const double& phi,
                          const arma::mat& cor_matrix,
                          const arma::vec& weights_vector) {
  const arma::vec variance_vector = variance(family, mu_vector);
  const arma::vec sd_vector = arma::sqrt(variance_vector / weights_vector);
  const arma::uword m = sd_vector.n_elem;

  if (m == 1) {
    arma::mat ans(1, 1);
    ans(0, 0) = phi * sd_vector[0] * sd_vector[0];
    return ans;
  }

  arma::mat ans;
  if (is_contiguous_1based(repeated_vector)) {
    const arma::uword r0 = static_cast<arma::uword>(repeated_vector[0]) - 1;
    const arma::uword r1 = r0 + m - 1;
    ans = cor_matrix.submat(r0, r0, r1, r1);
  } else {
    ans = subset_matrix(cor_matrix, repeated_vector);
  }

  ans.each_col() %= sd_vector;
  ans.each_row() %= sd_vector.t();
  ans *= phi;

  return ans;
}
//==============================================================================
