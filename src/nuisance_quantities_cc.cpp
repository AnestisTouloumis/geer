#include <RcppArmadillo.h>
#include "variance_functions.h"
#include "utils.h"
using namespace Rcpp;


//============================ pearson residuals ===============================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec get_pearson_residuals(const char * family,
                                const arma::vec & y_vector,
                                const arma::vec & mu_vector,
                                const arma::vec & weights_vector) {
  arma::vec ans =
    (y_vector - mu_vector) % sqrt(weights_vector/variance(family, mu_vector));
  return(ans);
}
//==============================================================================


//============================ phi hat =========================================
// [[Rcpp::export]]
double get_phi_hat(const arma::vec & pearson_residuals_vector,
                   const int & params_no) {
  int obs_no = pearson_residuals_vector.n_elem;
  double ans = arma::accu(pow(pearson_residuals_vector, 2))/(obs_no - params_no);
  if(ans < DBL_EPSILON) ans = 10 * DBL_EPSILON;
  return(ans);
}
//==============================================================================


//============================ exchangeable alpha hat ==========================
// [[Rcpp::export]]
double alpha_hat_exchangeable(const arma::vec & pearson_residuals_vector,
                              const arma::vec & id_vector,
                              const double & phi,
                              const int & params_no) {
  int sample_size = max(id_vector);
  double num = 0;
  double den = 0;
  for(int i = 1; i < sample_size + 1; i++) {
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec pearson_residuals_vector_i = pearson_residuals_vector(id_vector_i);
    int cluster_size_i = id_vector_i.n_elem;
    for(int j = 1; j < cluster_size_i; j++) {
      for(int k = j + 1; k < cluster_size_i + 1; k++) {
        num +=
          pearson_residuals_vector_i(j - 1) * pearson_residuals_vector_i(k - 1);
      }
    }
    den += cluster_size_i * (cluster_size_i - 1) * 0.5;
  }
  double ans = num/((den - params_no) * phi);
  return(ans);
}
//==============================================================================


//============================ ar1 alpha hat ===================================
// [[Rcpp::export]]
double alpha_hat_ar1(const arma::vec & pearson_residuals_vector,
                     const arma::vec & id_vector,
                     const arma::vec & repeated_vector,
                     const double & phi,
                     const int & params_no) {
  int sample_size = max(id_vector);
  double num = 0;
  double den = 0;
  for(int i = 1; i < sample_size + 1; i++) {
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec pearson_residuals_vector_i = pearson_residuals_vector(id_vector_i);
    arma::vec time_diff_vector_i = arma::diff(repeated_vector(id_vector_i));
    int cluster_size_i = id_vector_i.n_elem;
    for(int j = 1; j < cluster_size_i; j++) {
      if(time_diff_vector_i(j - 1) == 1) {
        num +=
          pearson_residuals_vector_i(j - 1) * pearson_residuals_vector_i(j);
        den += 1;
      }
    }
  }
  double ans = num/((den - params_no) * phi);
  return(ans);
}
//==============================================================================


//============================ unstructured alpha hat ==========================
// [[Rcpp::export]]
arma::vec alpha_hat_unstructured(const arma::vec & pearson_residuals_vector,
                                 const arma::vec & id_vector,
                                 const arma::vec & repeated_vector,
                                 const double & phi,
                                 const int & params_no) {
  int sample_size = max(id_vector);
  int time_max = max(repeated_vector);
  double time_pairs = time_max * (time_max - 1) / 2;
  arma::vec num(time_pairs);
  arma::vec den(time_pairs);
  for(int i = 1; i < sample_size + 1; i++) {
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec repeated_vector_i = repeated_vector(id_vector_i);
    arma::vec pearson_residuals_vector_i = pearson_residuals_vector(id_vector_i);
    int cluster_size_i = id_vector_i.n_elem;
    for(int j = 1; j < cluster_size_i; j++) {
      int index_j = repeated_vector_i(j - 1);
      double combn_j = index_j * (index_j + 1) / 2;
      for(int k = j + 1; k < cluster_size_i + 1; k++) {
        int time_index = time_max * (index_j - 1) + repeated_vector_i(k - 1) - combn_j;
        num(time_index - 1) +=
          pearson_residuals_vector_i(j - 1) * pearson_residuals_vector_i(k - 1);
        den(time_index - 1) += 1;
      }
    }
  }
  arma::vec ans = num/((den - params_no) * phi);
  return(ans);
}
//==============================================================================


//============================ m-dependent alpha hat ===========================
// [[Rcpp::export]]
arma::vec alpha_hat_mdependent(const arma::vec & pearson_residuals_vector,
                               const arma::vec & id_vector,
                               const arma::vec & repeated_vector,
                               const double & phi,
                               const int & params_no,
                               const int & mdependence) {
  int sample_size = max(id_vector);
  arma::vec num(mdependence);
  arma::vec den(mdependence);
  for(int i = 1; i < sample_size + 1; i++) {
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec repeated_vector_i = repeated_vector(id_vector_i);
    arma::vec pearson_residuals_vector_i = pearson_residuals_vector(id_vector_i);
    int cluster_size_i = id_vector_i.n_elem;
    for(int j = 1; j < cluster_size_i; j++) {
      int index_j = repeated_vector_i(j - 1);
      for(int k = j + 1; k < cluster_size_i + 1; k++) {
        int index_k = repeated_vector_i(k - 1);
        int diff_int = index_k - index_j;
        if(diff_int < mdependence + 1) {
          num(diff_int - 1) +=
            pearson_residuals_vector_i(j - 1) * pearson_residuals_vector_i(k - 1);
          den(diff_int - 1) += 1;
        }
      }
    }
  }
  arma::vec ans = num/((den - params_no) * phi);
  return(ans);
}
//==============================================================================


//============================ toeplitz alpha hat ==============================
// [[Rcpp::export]]
arma::vec alpha_hat_toeplitz(const arma::vec & pearson_residuals_vector,
                             const arma::vec & id_vector,
                             const arma::vec & repeated_vector,
                             const double & phi,
                             const int & params_no) {
  int sample_size = max(id_vector);
  int time_max = max(repeated_vector);
  arma::vec num(time_max - 1);
  arma::vec den(time_max - 1);
  for(int i = 1; i < sample_size + 1; i++) {
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec repeated_vector_i = repeated_vector(id_vector_i);
    arma::vec pearson_residuals_vector_i = pearson_residuals_vector(id_vector_i);
    int cluster_size_i = id_vector_i.n_elem;
    for(int j = 1; j < cluster_size_i; j++) {
      int index_j = repeated_vector_i(j - 1);
      for(int k = j + 1; k < cluster_size_i + 1; k++) {
        int index_k = repeated_vector_i(k - 1);
        num(index_k - index_j - 1) +=
          pearson_residuals_vector_i(j - 1) * pearson_residuals_vector_i(k - 1);
        den(index_k - index_j - 1) += 1;
      }
    }
  }
  arma::vec ans = num/((den - params_no) * phi);
  return(ans);
}
//==============================================================================


//============================ alpha hat =======================================
// [[Rcpp::export]]
arma::vec get_alpha_hat(const char * correlation_structure,
                        const arma::vec & pearson_residuals_vector,
                        const arma::vec & id_vector,
                        const arma::vec & repeated_vector,
                        const double & phi,
                        const int & params_no,
                        const int & mdependence) {
  arma::vec ans;
  if(std::strcmp(correlation_structure, "independence") == 0){
    ans.fill(0.0);
  }else if(std::strcmp(correlation_structure, "exchangeable") == 0){
    ans = alpha_hat_exchangeable(pearson_residuals_vector,
                                 id_vector,
                                 phi,
                                 params_no);
  }else if(std::strcmp(correlation_structure, "ar1") == 0){
    ans = alpha_hat_ar1(pearson_residuals_vector,
                        id_vector,
                        repeated_vector,
                        phi,
                        params_no);
  }else if(std::strcmp(correlation_structure, "m-dependent") == 0){
    ans = alpha_hat_mdependent(pearson_residuals_vector,
                               id_vector,
                               repeated_vector,
                               phi,
                               params_no,
                               mdependence);
  }else if(std::strcmp(correlation_structure, "unstructured") == 0){
    ans = alpha_hat_unstructured(pearson_residuals_vector,
                                 id_vector,
                                 repeated_vector,
                                 phi,
                                 params_no);
  }else if(std::strcmp(correlation_structure, "toeplitz") == 0){
    ans = alpha_hat_toeplitz(pearson_residuals_vector,
                             id_vector,
                             repeated_vector,
                             phi,
                             params_no);
  }
  double upper_bound = 1 - 10 * DBL_EPSILON;
  double lower_bound = 10 * DBL_EPSILON - 1;
  int ans_length = ans.n_elem;
  for(int i = 0; i < ans_length; i++){
    if(ans[i] >= upper_bound) ans[i] = upper_bound;
    if(ans[i] <= lower_bound) ans[i] = lower_bound;
  }
  return(ans);
}
//==============================================================================


//============================ independence ====================================
// [[Rcpp::export]]
arma::mat correlation_independence(const int & dimension) {
  arma::mat ans = arma::eye(dimension, dimension);
  return(ans);
}
//==============================================================================


//============================ exchangeable ====================================
// [[Rcpp::export]]
arma::mat correlation_exchangeable(const arma::vec & alpha_vector,
                                   const int & dimension) {
  arma::mat ans(dimension, dimension);
  ans.fill(alpha_vector(0));
  ans.diag().fill(1);
  return(ans);
}
//==============================================================================


//============================ ar1 =============================================
// [[Rcpp::export]]
arma::mat correlation_ar1(const arma::vec & alpha_vector,
                          const int & dimension) {
  arma::vec ans_elements(dimension);
  for(int i = 0; i < dimension; ++i) {
    ans_elements(i) = std::pow(alpha_vector(0), i);
  }
  arma::mat ans = arma::toeplitz(ans_elements);
  return(ans);
}
//==============================================================================


//============================ m-dependent =====================================
// [[Rcpp::export]]
arma::mat correlation_mdependent(const arma::vec & alpha_vector,
                                 const int & dimension) {
  arma::vec ones_vector(1);
  ones_vector.ones();
  int k = alpha_vector.n_elem;
  arma::vec zeros_vector(dimension - k - 1);
  arma::vec toeplitz_vector = arma::join_cols(ones_vector,
                                              alpha_vector,
                                              zeros_vector);
  arma::mat ans = arma::toeplitz(toeplitz_vector);
  return(ans);
}
//==============================================================================


//============================ toeplitz ========================================
// [[Rcpp::export]]
arma::mat correlation_toeplitz(const arma::vec & alpha_vector) {
  arma::vec ones_vector(1);
  ones_vector.ones();
  arma::vec toeplitz_vector = arma::join_cols(ones_vector, alpha_vector);
  arma::mat ans = arma::toeplitz(toeplitz_vector);
  return(ans);
}
//==============================================================================


//============================ unstructured ====================================
// [[Rcpp::export]]
arma::mat correlation_unstructured(const arma::vec & alpha_vector,
                                   const int & dimension) {
  arma::mat ans_lt = arma::eye(dimension, dimension);
  arma::uvec lt_indices = arma::trimatl_ind(arma::size(ans_lt), - 1);
  ans_lt.elem(lt_indices) = alpha_vector;
  arma::mat ans = symmatl(ans_lt);
  return(ans);
}
//==============================================================================


//============================ correlation matrix given rho vector =============
// [[Rcpp::export]]
arma::mat get_correlation_matrix(const char * correlation_structure,
                                 const arma::vec & alpha_vector,
                                 const int & dimension) {
  arma::mat ans(dimension, dimension);
  if(std::strcmp(correlation_structure, "independence") == 0){
    ans = correlation_independence(dimension);
  }else if(std::strcmp(correlation_structure, "exchangeable") == 0){
    ans = correlation_exchangeable(alpha_vector, dimension);
  }else if(std::strcmp(correlation_structure, "ar1") == 0){
    ans = correlation_ar1(alpha_vector, dimension);
  }else if(std::strcmp(correlation_structure, "m-dependent") == 0){
    ans = correlation_mdependent(alpha_vector, dimension);
  }else if(std::strcmp(correlation_structure, "unstructured") == 0){
    ans = correlation_unstructured(alpha_vector, dimension);
  }else if(std::strcmp(correlation_structure, "toeplitz") == 0){
    ans = correlation_toeplitz(alpha_vector);
  }else if(std::strcmp(correlation_structure, "fixed") == 0){
    ans = correlation_unstructured(alpha_vector, dimension);
  }
  return(ans);
}
//==============================================================================


// =========================== subject-specific weight matrix ==================
// [[Rcpp::export]]
arma::mat get_v_matrix_cc(const char * family,
                          const arma::vec & mu_vector,
                          const arma::vec & repeated_vector,
                          const double & phi,
                          const arma::mat & cor_matrix,
                          const arma::vec & weights_vector) {
  arma::vec sd_vector =
    sqrt(variance(family, mu_vector)/weights_vector);
  if(repeated_vector.n_elem == 1) {
    return(phi * (sd_vector * trans(sd_vector)));
  } else {
    arma::mat ans =
      phi *
      subset_matrix(cor_matrix, repeated_vector) %
      (sd_vector * trans(sd_vector));
    return(ans);
  }
}
//==============================================================================
