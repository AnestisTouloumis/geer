#include "nuisance_estimation_cc.h"
#include "utils.h"
#include "variance_functions.h"
#include <cstring>
#include <vector>


//============================ pearson residuals (enum) ========================
arma::vec get_pearson_residuals(FamilyCode family_code,
                                const arma::vec& y_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& weights_vector) {
  const arma::vec scale = arma::sqrt(weights_vector / variance_function(family_code, mu_vector));
  return (y_vector - mu_vector) % scale;
}
//==============================================================================


//============================ pearson residuals (char*) =======================
// [[Rcpp::export]]
arma::vec get_pearson_residuals(const char* family,
                                const arma::vec& y_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& weights_vector) {
  reset_lu_fallback_warning();
  return get_pearson_residuals(parse_family(family), y_vector, mu_vector, weights_vector);
}
//==============================================================================


//============================ phi hat =========================================
double get_phi_hat_cc(const arma::vec& pearson_residuals_vector,
                   const int params_no) {
  const double n = static_cast<double>(pearson_residuals_vector.n_elem);
  const double denominator = n - static_cast<double>(params_no);
  if (denominator <= 0.0) {
    Rcpp::stop("get_phi_hat_cc failed: the dispersion cannot be "
      "estimated: too few observations relative to the "
      "number "
               "of regression parameters");
  }
  const double result =
    arma::accu(arma::square(pearson_residuals_vector)) / denominator;
  if (!R_FINITE(result) || result <= 0.0) {
    Rcpp::stop(
      "the dispersion estimate (%.6g) is not a usable "
      "positive value; the fitted means may "
      "reproduce the data exactly",
      result
    );
  }
  return result;
}
//==============================================================================


namespace {

//============================ exchangeable alpha hat ==========================
double get_alpha_hat_exchangeable(const arma::vec& pearson_residuals_vector,
                              const std::vector<Cluster>& clusters,
                              const double phi,
                              const int params_no) {
  double num = 0.0;
  double den = 0.0;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.first_row;
      const arma::uword last_row = cluster.after_last_row - 1;
      const arma::vec pearson_residuals_vector_i =
        pearson_residuals_vector.subvec(first_row, last_row);
      const arma::uword cluster_size_i =
        cluster.after_last_row - cluster.first_row;
      for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
        for (arma::uword k = j + 1; k < cluster_size_i; ++k) {
          num += pearson_residuals_vector_i[j] * pearson_residuals_vector_i[k];
        }
      }
      den += static_cast<double>(cluster_size_i) *
             static_cast<double>(cluster_size_i - 1) * 0.5;
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_alpha_hat_exchangeable", cluster_index,
                                   cluster.id_value, e);
    }
  }
  const double denominator = (den - static_cast<double>(params_no)) * phi;
  if (denominator <= 0.0) {
    Rcpp::stop(
      "the exchangeable correlation cannot be estimated: "
      "too few pairs of observations relative "
      "to the number of regression parameters"
    );
  }
  return num / denominator;
}
//==============================================================================


//============================ ar1 alpha hat ===================================
double get_alpha_hat_ar1(const arma::vec& pearson_residuals_vector,
                     const std::vector<Cluster>& clusters,
                     const arma::vec& repeated_vector,
                     const double phi,
                     const int params_no) {
  double num = 0.0;
  double den = 0.0;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.first_row;
      const arma::uword last_row = cluster.after_last_row - 1;
      const arma::vec pearson_residuals_vector_i =
        pearson_residuals_vector.subvec(first_row, last_row);
      const arma::vec repeated_vector_i = repeated_vector.subvec(first_row, last_row);
      const arma::uword cluster_size_i =
        cluster.after_last_row - cluster.first_row;
      for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
        if (repeated_vector_i[j + 1] - repeated_vector_i[j] == 1.0) {
          num += pearson_residuals_vector_i[j] * pearson_residuals_vector_i[j + 1];
          den += 1.0;
        }
      }
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_alpha_hat_ar1", cluster_index,
                                   cluster.id_value, e);
    }
  }
  const double denominator_ar1 = (den - static_cast<double>(params_no)) * phi;
  if (denominator_ar1 <= 0.0) {
    Rcpp::stop(
      "the AR(1) correlation cannot be estimated: too few "
      "pairs of consecutive observations "
      "relative to the number of regression parameters"
    );
  }
  return num / denominator_ar1;
}
//==============================================================================


//============================ unstructured alpha hat ==========================
arma::vec get_alpha_hat_unstructured(const arma::vec& pearson_residuals_vector,
                                 const std::vector<Cluster>& clusters,
                                 const arma::vec& repeated_vector,
                                 const double phi,
                                 const int params_no) {
  const arma::uword time_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::uword time_pairs = time_max * (time_max - 1) / 2;
  arma::vec num(time_pairs, arma::fill::zeros);
  arma::vec den(time_pairs, arma::fill::zeros);
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.first_row;
      const arma::uword last_row = cluster.after_last_row - 1;
      const arma::vec repeated_vector_i = repeated_vector.subvec(first_row, last_row);
      const arma::vec pearson_residuals_vector_i =
        pearson_residuals_vector.subvec(first_row, last_row);
      const arma::uword cluster_size_i =
        cluster.after_last_row - cluster.first_row;
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
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_alpha_hat_unstructured", cluster_index,
                                   cluster.id_value, e);
    }
  }
  const arma::vec denominator_unstr = (den - static_cast<double>(params_no)) * phi;
  if (arma::any(denominator_unstr <= 0.0)) {
    Rcpp::stop(
      "the unstructured correlation cannot be estimated: "
      "some pairs of occasions have too few "
      "observations relative to the number of regression parameters"
    );
  }
  return num / denominator_unstr;
}
//==============================================================================


//============================ m-dependent alpha hat ===========================
arma::vec get_alpha_hat_mdependent(const arma::vec& pearson_residuals_vector,
                               const std::vector<Cluster>& clusters,
                               const arma::vec& repeated_vector,
                               const double phi,
                               const int params_no,
                               const int mdependence) {
  arma::vec num(mdependence, arma::fill::zeros);
  arma::vec den(mdependence, arma::fill::zeros);
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.first_row;
      const arma::uword last_row = cluster.after_last_row - 1;
      const arma::vec repeated_vector_i = repeated_vector.subvec(first_row, last_row);
      const arma::vec pearson_residuals_vector_i =
        pearson_residuals_vector.subvec(first_row, last_row);
      const arma::uword cluster_size_i =
        cluster.after_last_row - cluster.first_row;
      for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
        const arma::uword index_j =
          static_cast<arma::uword>(repeated_vector_i[j]);
        for (arma::uword k = j + 1; k < cluster_size_i; ++k) {
          const arma::uword index_k =
            static_cast<arma::uword>(repeated_vector_i[k]);
          if (index_k <= index_j) {
            signal_input_error(
      "'repeated' must be strictly increasing within each "
      "cluster");
          }
          const arma::uword diff_int = index_k - index_j;

          if (diff_int < static_cast<arma::uword>(mdependence + 1)) {
            num[diff_int - 1] +=
              pearson_residuals_vector_i[j] * pearson_residuals_vector_i[k];
            den[diff_int - 1] += 1.0;
          }
        }
      }
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_alpha_hat_mdependent", cluster_index,
                                   cluster.id_value, e);
    }
  }
  const arma::vec denominator_mdep = (den - static_cast<double>(params_no)) * phi;
  if (arma::any(denominator_mdep <= 0.0)) {
    Rcpp::stop(
      "the m-dependent correlation cannot be estimated: "
      "some lags have too few observations "
      "relative to the number of regression parameters"
    );
  }
  return num / denominator_mdep;
}
//==============================================================================


//============================ toeplitz alpha hat ==============================
arma::vec get_alpha_hat_toeplitz(const arma::vec& pearson_residuals_vector,
                             const std::vector<Cluster>& clusters,
                             const arma::vec& repeated_vector,
                             const double phi,
                             const int params_no) {
  const arma::uword time_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  if (time_max < 2) {
    Rcpp::stop(
      "the Toeplitz correlation cannot be estimated: no "
      "cluster contains a pair of observations"
    );
  }
  arma::vec num(time_max - 1, arma::fill::zeros);
  arma::vec den(time_max - 1, arma::fill::zeros);
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.first_row;
      const arma::uword last_row = cluster.after_last_row - 1;
      const arma::vec repeated_vector_i = repeated_vector.subvec(first_row, last_row);
      const arma::vec pearson_residuals_vector_i =
        pearson_residuals_vector.subvec(first_row, last_row);
      const arma::uword cluster_size_i =
        cluster.after_last_row - cluster.first_row;
      for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
        const arma::uword index_j =
          static_cast<arma::uword>(repeated_vector_i[j]);
        for (arma::uword k = j + 1; k < cluster_size_i; ++k) {
          const arma::uword index_k =
            static_cast<arma::uword>(repeated_vector_i[k]);
          const arma::uword lag = index_k - index_j;
          if (index_k <= index_j || lag > num.n_elem) {
            signal_input_error(
      "'repeated' must be strictly increasing within each "
      "cluster");
          }
          num[lag - 1] +=
            pearson_residuals_vector_i[j] * pearson_residuals_vector_i[k];
          den[lag - 1] += 1.0;
        }
      }
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("get_alpha_hat_toeplitz", cluster_index,
                                   cluster.id_value, e);
    }
  }
  const arma::vec denominator_toep = (den - static_cast<double>(params_no)) * phi;
  if (arma::any(denominator_toep <= 0.0)) {
    Rcpp::stop(
      "the Toeplitz correlation cannot be estimated: some "
      "lags have too few observations "
      "relative to the number of regression parameters"
    );
  }
  return num / denominator_toep;
}
//==============================================================================


}  // namespace

//============================ alpha hat =======================================
arma::vec get_alpha_hat_cc(const char* correlation_structure,
                        const arma::vec& pearson_residuals_vector,
                        const std::vector<Cluster>& clusters,
                        const arma::vec& repeated_vector,
                        const double phi,
                        const int params_no,
                        const int mdependence) {
  if (std::strcmp(correlation_structure, "independence") == 0) {
    return arma::vec();
  }
  arma::vec result;
  if (std::strcmp(correlation_structure, "exchangeable") == 0) {
    result.set_size(1);
    result[0] = get_alpha_hat_exchangeable(pearson_residuals_vector,
                                       clusters,
                                       phi,
                                       params_no);
  } else if (std::strcmp(correlation_structure, "ar1") == 0) {
    result.set_size(1);
    result[0] = get_alpha_hat_ar1(pearson_residuals_vector,
                              clusters,
                              repeated_vector,
                              phi,
                              params_no);
  } else if (std::strcmp(correlation_structure, "m-dependent") == 0) {
    result = get_alpha_hat_mdependent(pearson_residuals_vector,
                                  clusters,
                                  repeated_vector,
                                  phi,
                                  params_no,
                                  mdependence);
  } else if (std::strcmp(correlation_structure, "unstructured") == 0) {
    result = get_alpha_hat_unstructured(pearson_residuals_vector,
                                    clusters,
                                    repeated_vector,
                                    phi,
                                    params_no);
  } else if (std::strcmp(correlation_structure, "toeplitz") == 0) {
    result = get_alpha_hat_toeplitz(pearson_residuals_vector,
                                clusters,
                                repeated_vector,
                                phi,
                                params_no);
  } else if (std::strcmp(correlation_structure, "fixed") == 0) {
    signal_input_error(
      "'alpha_fixed' must be 1 for the \"fixed\" correlation "
      "structure, so that "
      "alpha is not re-estimated"
    );
  } else {
    signal_input_error("'correlation_structure' must be one of: "
      "independence, exchangeable, ar1, m-dependent, "
      "unstructured, toeplitz, fixed; received '%s'",
               correlation_structure);
  }
  return result;
}
//==============================================================================
