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
  const arma::vec scale = arma::sqrt(weights_vector / geer::variance_function(family_code, mu_vector));
  return (y_vector - mu_vector) % scale;
}
//==============================================================================


//============================ pearson residuals (char*) =======================
// [[Rcpp::export]]
arma::vec get_pearson_residuals(const char* family,
                                const arma::vec& y_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& weights_vector) {
  return get_pearson_residuals(parse_family(family), y_vector, mu_vector, weights_vector);
}
//==============================================================================


//============================ phi hat =========================================
double get_phi_hat(const arma::vec& pearson_residuals_vector,
                   const int params_no) {
  const double n = static_cast<double>(pearson_residuals_vector.n_elem);
  const double denominator = n - static_cast<double>(params_no);
  if (denominator <= 0.0) {
    Rcpp::stop("get_phi_hat: non-positive denominator.");
  }
  const double result =
    arma::accu(arma::square(pearson_residuals_vector)) / denominator;
  // A dispersion estimate that is zero (or not finite) means the fitted means
  // reproduce the data, so the scaled residuals and the working covariance
  // matrices are undefined. Reporting it is preferred to flooring the value,
  // which would silently inflate every scaled quantity that divides by it. The
  // test is deliberately not an absolute threshold: phi scales with the square
  // of the response, so small but valid values occur for responses measured in
  // small units.
  if (!R_FINITE(result) || result <= 0.0) {
    Rcpp::stop(
      "get_phi_hat: the dispersion estimate is not a usable positive value "
      "(%.6g); the fitted means may reproduce the data exactly.",
      result
    );
  }
  return result;
}
//==============================================================================


namespace {

//============================ exchangeable alpha hat ==========================
double alpha_hat_exchangeable(const arma::vec& pearson_residuals_vector,
                              const std::vector<Cluster>& clusters,
                              const double phi,
                              const int params_no) {
  double num = 0.0;
  double den = 0.0;
  for (const auto& cluster : clusters) {
    const arma::uword first_row = cluster.start;
    const arma::uword last_row = cluster.end - 1;
    const arma::vec pearson_residuals_vector_i =
      pearson_residuals_vector.subvec(first_row, last_row);
    const arma::uword cluster_size_i = cluster.end - cluster.start;
    for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
      for (arma::uword k = j + 1; k < cluster_size_i; ++k) {
        num += pearson_residuals_vector_i[j] * pearson_residuals_vector_i[k];
      }
    }
    den += static_cast<double>(cluster_size_i) *
      static_cast<double>(cluster_size_i - 1) * 0.5;
  }
  const double denominator = (den - static_cast<double>(params_no)) * phi;
  if (denominator <= 0.0) {
    Rcpp::stop(
      "alpha_hat_exchangeable: non-positive denominator -- "
      "too few observation pairs relative to the number of parameters."
    );
  }
  return num / denominator;
}
//==============================================================================


//============================ ar1 alpha hat ===================================
double alpha_hat_ar1(const arma::vec& pearson_residuals_vector,
                     const std::vector<Cluster>& clusters,
                     const arma::vec& repeated_vector,
                     const double phi,
                     const int params_no) {
  double num = 0.0;
  double den = 0.0;
  for (const auto& cluster : clusters) {
    const arma::uword first_row = cluster.start;
    const arma::uword last_row = cluster.end - 1;
    const arma::vec pearson_residuals_vector_i =
      pearson_residuals_vector.subvec(first_row, last_row);
    const arma::vec repeated_vector_i = repeated_vector.subvec(first_row, last_row);
    const arma::uword cluster_size_i = cluster.end - cluster.start;
    for (arma::uword j = 0; j + 1 < cluster_size_i; ++j) {
      if (repeated_vector_i[j + 1] - repeated_vector_i[j] == 1.0) {
        num += pearson_residuals_vector_i[j] * pearson_residuals_vector_i[j + 1];
        den += 1.0;
      }
    }
  }
  const double denominator_ar1 = (den - static_cast<double>(params_no)) * phi;
  if (denominator_ar1 <= 0.0) {
    Rcpp::stop(
      "alpha_hat_ar1: non-positive denominator -- "
      "too few consecutive observation pairs relative to the number of parameters."
    );
  }
  return num / denominator_ar1;
}
//==============================================================================


//============================ unstructured alpha hat ==========================
arma::vec alpha_hat_unstructured(const arma::vec& pearson_residuals_vector,
                                 const std::vector<Cluster>& clusters,
                                 const arma::vec& repeated_vector,
                                 const double phi,
                                 const int params_no) {
  const arma::uword time_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::uword time_pairs = time_max * (time_max - 1) / 2;
  arma::vec num(time_pairs, arma::fill::zeros);
  arma::vec den(time_pairs, arma::fill::zeros);
  for (const auto& cluster : clusters) {
    const arma::uword first_row = cluster.start;
    const arma::uword last_row = cluster.end - 1;
    const arma::vec repeated_vector_i = repeated_vector.subvec(first_row, last_row);
    const arma::vec pearson_residuals_vector_i =
      pearson_residuals_vector.subvec(first_row, last_row);
    const arma::uword cluster_size_i = cluster.end - cluster.start;
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
  const arma::vec denominator_unstr = (den - static_cast<double>(params_no)) * phi;
  if (arma::any(denominator_unstr <= 0.0)) {
    Rcpp::stop(
      "alpha_hat_unstructured: non-positive denominator for one or more time pairs -- "
      "too few observations relative to the number of parameters."
    );
  }
  return num / denominator_unstr;
}
//==============================================================================


//============================ m-dependent alpha hat ===========================
arma::vec alpha_hat_mdependent(const arma::vec& pearson_residuals_vector,
                               const std::vector<Cluster>& clusters,
                               const arma::vec& repeated_vector,
                               const double phi,
                               const int params_no,
                               const int mdependence) {
  arma::vec num(mdependence, arma::fill::zeros);
  arma::vec den(mdependence, arma::fill::zeros);
  for (const auto& cluster : clusters) {
    const arma::uword first_row = cluster.start;
    const arma::uword last_row = cluster.end - 1;
    const arma::vec repeated_vector_i = repeated_vector.subvec(first_row, last_row);
    const arma::vec pearson_residuals_vector_i =
      pearson_residuals_vector.subvec(first_row, last_row);
    const arma::uword cluster_size_i = cluster.end - cluster.start;
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
  const arma::vec denominator_mdep = (den - static_cast<double>(params_no)) * phi;
  if (arma::any(denominator_mdep <= 0.0)) {
    Rcpp::stop(
      "alpha_hat_mdependent: non-positive denominator for one or more lags -- "
      "too few observations relative to the number of parameters."
    );
  }
  return num / denominator_mdep;
}
//==============================================================================


//============================ toeplitz alpha hat ==============================
arma::vec alpha_hat_toeplitz(const arma::vec& pearson_residuals_vector,
                             const std::vector<Cluster>& clusters,
                             const arma::vec& repeated_vector,
                             const double phi,
                             const int params_no) {
  const arma::uword time_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec num(time_max - 1, arma::fill::zeros);
  arma::vec den(time_max - 1, arma::fill::zeros);
  for (const auto& cluster : clusters) {
    const arma::uword first_row = cluster.start;
    const arma::uword last_row = cluster.end - 1;
    const arma::vec repeated_vector_i = repeated_vector.subvec(first_row, last_row);
    const arma::vec pearson_residuals_vector_i =
      pearson_residuals_vector.subvec(first_row, last_row);
    const arma::uword cluster_size_i = cluster.end - cluster.start;
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
  const arma::vec denominator_toep = (den - static_cast<double>(params_no)) * phi;
  if (arma::any(denominator_toep <= 0.0)) {
    Rcpp::stop(
      "alpha_hat_toeplitz: non-positive denominator for one or more lags -- "
      "too few observations relative to the number of parameters."
    );
  }
  return num / denominator_toep;
}
//==============================================================================


} // namespace

//============================ alpha hat =======================================
arma::vec get_alpha_hat(const char* correlation_structure,
                        const arma::vec& pearson_residuals_vector,
                        const std::vector<Cluster>& clusters,
                        const arma::vec& repeated_vector,
                        const double phi,
                        const int params_no,
                        const int mdependence) {
  if (std::strcmp(correlation_structure, "independence") == 0) {
    return arma::vec(); // no alpha parameters for independence
  }
  arma::vec result;
  if (std::strcmp(correlation_structure, "exchangeable") == 0) {
    result.set_size(1);
    result[0] = alpha_hat_exchangeable(pearson_residuals_vector,
                                    clusters,
                                    phi,
                                    params_no);
  } else if (std::strcmp(correlation_structure, "ar1") == 0) {
    result.set_size(1);
    result[0] = alpha_hat_ar1(pearson_residuals_vector,
                           clusters,
                           repeated_vector,
                           phi,
                           params_no);
  } else if (std::strcmp(correlation_structure, "m-dependent") == 0) {
    result = alpha_hat_mdependent(pearson_residuals_vector,
                               clusters,
                               repeated_vector,
                               phi,
                               params_no,
                               mdependence);
  } else if (std::strcmp(correlation_structure, "unstructured") == 0) {
    result = alpha_hat_unstructured(pearson_residuals_vector,
                                 clusters,
                                 repeated_vector,
                                 phi,
                                 params_no);
  } else if (std::strcmp(correlation_structure, "toeplitz") == 0) {
    result = alpha_hat_toeplitz(pearson_residuals_vector,
                             clusters,
                             repeated_vector,
                             phi,
                             params_no);
  } else if (std::strcmp(correlation_structure, "fixed") == 0) {
    Rcpp::stop(
      "get_alpha_hat: correlation structure is \"fixed\" but alpha_fixed must "
      "be 1 when using a fixed correlation -- alpha should not be re-estimated."
    );
  } else {
    Rcpp::stop("get_alpha_hat: unsupported correlation structure \"%s\".",
               correlation_structure);
  }
  return result;
}
//==============================================================================
