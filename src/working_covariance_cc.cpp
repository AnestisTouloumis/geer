#include "working_covariance_cc.h"
#include "utils.h"
#include "variance_functions.h"
#include <cmath>
#include <cstring>

namespace {
bool is_contiguous_1based(const arma::vec& r) {
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


//============================ independence ====================================
arma::mat get_correlation_matrix_independence(const arma::uword dimension) {
  return arma::eye(dimension, dimension);
}
//==============================================================================


//============================ exchangeable ====================================
arma::mat get_correlation_matrix_exchangeable(const arma::vec& alpha_vector,
                                   const arma::uword dimension) {
  arma::mat result(dimension, dimension);
  result.fill(alpha_vector[0]);
  result.diag().fill(1.0);
  return result;
}
//==============================================================================


//============================ ar1 =============================================
arma::mat get_correlation_matrix_ar1(const arma::vec& alpha_vector,
                          const arma::uword dimension) {
  arma::vec result(dimension, arma::fill::zeros);
  result[0] = 1.0;

  for (arma::uword i = 1; i < dimension; ++i) {
    result[i] = result[i - 1] * alpha_vector[0];
  }
  return arma::toeplitz(result);
}
//==============================================================================


//============================ m-dependent =====================================
arma::mat get_correlation_matrix_mdependent(const arma::vec& alpha_vector,
                                 const arma::uword dimension) {
  const arma::uword k = alpha_vector.n_elem;

  if (k + 1 > dimension) {
    Rcpp::stop("get_correlation_matrix_mdependent: 'alpha_vector' is too long");
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
arma::mat get_correlation_matrix_toeplitz(const arma::vec& alpha_vector,
                               const arma::uword dimension) {
  if (!alpha_vector.is_empty() && alpha_vector.n_elem >= dimension) {
    Rcpp::stop("get_correlation_matrix_toeplitz: 'alpha_vector' is too long");
  }
  arma::vec toeplitz_vector(dimension, arma::fill::zeros);
  toeplitz_vector[0] = 1.0;
  if (!alpha_vector.is_empty()) {
    toeplitz_vector.subvec(1, alpha_vector.n_elem) = alpha_vector;
  }

  return arma::toeplitz(toeplitz_vector);
}
//==============================================================================


//============================ unstructured ====================================
arma::mat get_correlation_matrix_unstructured(const arma::vec& alpha_vector,
                                   const arma::uword dimension) {
  arma::mat ans_lt = arma::eye(dimension, dimension);
  if (dimension < 2) {
    return ans_lt;
  }
  const arma::uvec lt_indices = arma::trimatl_ind(arma::size(ans_lt), -1);
  ans_lt.elem(lt_indices) = alpha_vector;
  return arma::symmatl(ans_lt);
}
//==============================================================================


//============================ scalar alpha check ==============================
void check_alpha_scalar(const char* correlation_structure,
                        const arma::vec& alpha_vector) {
  if (alpha_vector.n_elem < 1) {
    Rcpp::stop("the %s correlation parameter must be supplied",
               correlation_structure);
  }
  if (!R_FINITE(alpha_vector[0])) {
    Rcpp::stop("the %s correlation parameter must be finite, not %.6g",
               correlation_structure, alpha_vector[0]);
  }
}
//==============================================================================

}  // namespace

//============================ correlation matrix given rho vector =============
// [[Rcpp::export]]
arma::mat get_correlation_matrix_cc(const char* correlation_structure,
                                 const arma::vec& alpha_vector,
                                 const arma::uword dimension) {
  reset_lu_fallback_warning();
  if (std::strcmp(correlation_structure, "independence") == 0) {
    return get_correlation_matrix_independence(dimension);
  }
  if (std::strcmp(correlation_structure, "ar1") == 0) {
    check_alpha_scalar(correlation_structure, alpha_vector);
    if (alpha_vector[0] <= -1.0 || alpha_vector[0] >= 1.0) {
      Rcpp::stop(
        "the ar1 correlation parameter must lie in (-1, 1), not %.6g",
        alpha_vector[0]
      );
    }
    return get_correlation_matrix_ar1(alpha_vector, dimension);
  }
  if (std::strcmp(correlation_structure, "exchangeable") == 0) {
    check_alpha_scalar(correlation_structure, alpha_vector);
    const double lower = -1.0 / static_cast<double>(dimension - 1);
    if (alpha_vector[0] <= lower || alpha_vector[0] >= 1.0) {
      Rcpp::stop(
        "the exchangeable correlation parameter must lie in (%.6g, 1) for the largest cluster "
        "size of %d, not %.6g",
        lower,
        static_cast<int>(dimension),
        alpha_vector[0]
      );
    }
    return get_correlation_matrix_exchangeable(alpha_vector, dimension);
  }
  for (arma::uword i = 0; i < alpha_vector.n_elem; ++i) {
    if (!R_FINITE(alpha_vector[i])) {
      Rcpp::stop("the %s correlation parameter alpha[%d] must be finite",
                 correlation_structure, static_cast<int>(i + 1));
    }
    if (alpha_vector[i] <= -1.0 || alpha_vector[i] >= 1.0) {
      Rcpp::stop(
        "the %s correlation parameter alpha[%d] must lie in (-1, 1), not %.6g",
        correlation_structure,
        static_cast<int>(i + 1),
        alpha_vector[i]
      );
    }
  }
  arma::mat cor_mat;
  if (std::strcmp(correlation_structure, "m-dependent") == 0) {
    cor_mat = get_correlation_matrix_mdependent(alpha_vector, dimension);
  } else if (std::strcmp(correlation_structure, "toeplitz") == 0) {
    cor_mat = get_correlation_matrix_toeplitz(alpha_vector, dimension);
  } else if (std::strcmp(correlation_structure, "unstructured") == 0 ||
             std::strcmp(correlation_structure, "fixed") == 0) {
    cor_mat = get_correlation_matrix_unstructured(alpha_vector, dimension);
  } else {
    Rcpp::stop("unsupported correlation structure: %s",
               correlation_structure);
  }
  arma::mat chol_result;
  if (!arma::chol(chol_result, cor_mat)) {
    Rcpp::stop(
      "the %s working correlation matrix is not positive definite; consider a simpler "
      "correlation structure or fixing alpha",
      correlation_structure
    );
  }
  return cor_mat;
}
//==============================================================================


//============================ working covariance matrix V_i (enum) ============
arma::mat get_v_matrix_cc(FamilyCode family_code,
                          const arma::vec& mu_vector,
                          const arma::vec& repeated_vector,
                          const double phi,
                          const arma::mat& cor_matrix,
                          const arma::vec& weights_vector) {
  const arma::vec sd_vector =
    arma::sqrt(variance_function(family_code, mu_vector) / weights_vector);
  const arma::uword cluster_size = sd_vector.n_elem;
  if (cluster_size == 1) {
    arma::mat result(1, 1);
    result(0, 0) = phi * sd_vector[0] * sd_vector[0];
    return result;
  }
  arma::mat result;
  if (is_contiguous_1based(repeated_vector)) {
    const arma::uword r0 = static_cast<arma::uword>(repeated_vector[0]) - 1;
    const arma::uword r1 = r0 + cluster_size - 1;
    result = cor_matrix.submat(r0, r0, r1, r1);
  } else {
    result = subset_matrix(cor_matrix, repeated_vector);
  }
  result.each_col() %= sd_vector;
  result.each_row() %= sd_vector.t();
  result *= phi;
  return result;
}
//==============================================================================


//============================ working covariance matrix V_i (char*) ===========
arma::mat get_v_matrix_cc(const char* family,
                          const arma::vec& mu_vector,
                          const arma::vec& repeated_vector,
                          const double phi,
                          const arma::mat& cor_matrix,
                          const arma::vec& weights_vector) {
  return get_v_matrix_cc(parse_family(family),
                         mu_vector,
                         repeated_vector,
                         phi,
                         cor_matrix,
                         weights_vector);
}
//==============================================================================
