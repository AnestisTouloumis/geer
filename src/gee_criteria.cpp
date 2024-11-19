#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link_functions.h"
#include "variance_functions.h"
#include "nuisance_quantities.h"
#include "nuisance_quantities_oddsratio.h"
#include "utils.h"
using namespace Rcpp;

//============================ naive matrix inverse - independence =============
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat get_naive_matrix_inverse_independence(const arma::vec & y_vector,
                                                const arma::mat & model_matrix,
                                                const arma::vec & id_vector,
                                                const char* link,
                                                const char* family,
                                                const arma::vec & mu_vector,
                                                const arma::vec & eta_vector,
                                                const double & phi) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  arma::mat ans = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec variance_vector = variance(family, mu_vector);
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::mat d_matrix_i = arma::diagmat(delta_vector(id_vector_i)) *
      model_matrix.rows(id_vector_i);
    arma::mat alpha_matrix_inverse =
      arma::diagmat(1/variance_vector(id_vector_i));
    ans += trans(d_matrix_i) * alpha_matrix_inverse * d_matrix_i;
  }
  return ans/phi;
}
//==============================================================================


//============================ get_sc_criteria =================================
// [[Rcpp::export()]]
Rcpp::List get_gee_criteria_sc_cw(const arma::vec & y_vector,
                                  const arma::vec & id_vector,
                                  const arma::vec & repeated_vector,
                                  const char * family,
                                  const arma::vec & mu_vector,
                                  const char * correlation_structure,
                                  const arma::vec & alpha_vector,
                                  const double & phi) {
  double sample_size = max(id_vector);
  arma::vec sc_criterion = arma::zeros(1, 1);
  double sum_log_det_working_covariance_matrices = 0;
  arma::vec s_vector = y_vector - mu_vector;
  arma::mat correlation_matrix_inverse_full =
    get_correlation_matrix_inverse(correlation_structure,
                                   alpha_vector,
                                   max(repeated_vector));
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec s_vector_i = s_vector(id_vector_i);
    arma::mat working_covariance_matrix_i =
      get_weight_matrix(family,
                        mu_vector(id_vector_i),
                        repeated_vector(id_vector_i),
                        phi,
                        correlation_matrix_inverse_full);
    sc_criterion += trans(s_vector_i) * working_covariance_matrix_i * s_vector_i;
    sum_log_det_working_covariance_matrices += log(det(working_covariance_matrix_i));
  }
  Rcpp::List ans;
  ans["sc"] = sc_criterion;
  ans["gp"] = 0.5 * (sc_criterion - sum_log_det_working_covariance_matrices);
  return ans;
}
//==============================================================================


//============================ covariance matrices with odds ratios ============
// [[Rcpp::export()]]
Rcpp::List get_gee_criteria_sc_cw_or(const arma::vec & y_vector,
                                     const arma::vec & id_vector,
                                     const arma::vec & repeated_vector,
                                     const arma::vec & mu_vector,
                                     const arma::vec & alpha_vector) {
  double sample_size = max(id_vector);
  arma::vec sc_criterion = arma::zeros(1, 1);
  double sum_log_det_working_covariance_matrices = 0;
  arma::vec s_vector = y_vector - mu_vector;
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec mu_vector_i = mu_vector(id_vector_i);
    arma::vec s_vector_i = s_vector(id_vector_i);
    arma::mat working_covariance_matrix_i =
      get_weight_matrix_inverse_or(mu_vector_i,
                                   get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                                                    max(repeated_vector),
                                                                    alpha_vector));
    sc_criterion += trans(s_vector_i) * working_covariance_matrix_i * s_vector_i;
    sum_log_det_working_covariance_matrices += log(det(working_covariance_matrix_i));
  }
  Rcpp::List ans;
  ans["sc"] = sc_criterion;
  ans["gp"] = 0.5 * (sc_criterion - sum_log_det_working_covariance_matrices);
  return ans;
}
