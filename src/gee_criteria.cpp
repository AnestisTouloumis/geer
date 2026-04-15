#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link_functions.h"
#include "variance_functions.h"
#include "nuisance_quantities_cc.h"
#include "nuisance_quantities_or.h"
#include "cluster_utils.h"
#include "utils.h"

//============================ naive matrix inverse - independence =============
// [[Rcpp::export]]
arma::mat get_naive_matrix_inverse_independence(const arma::vec& y_vector,
                                                const arma::mat& model_matrix,
                                                const arma::vec& id_vector,
                                                const char* link,
                                                const char* family,
                                                const arma::vec& mu_vector,
                                                const arma::vec& eta_vector,
                                                const double& phi,
                                                const arma::vec& weights_vector) {
  const arma::uword params_no = model_matrix.n_cols;
  const arma::vec delta_vector = mueta(link, eta_vector);
  const arma::vec variance_vector = variance(family, mu_vector);
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::mat ans(params_no, params_no, arma::fill::zeros);
  arma::mat d_matrix_i;
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::uword m = cl.end - cl.start;
    if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
      d_matrix_i.set_size(m, params_no);
    }
    d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector.subvec(a, b);
    const arma::vec scale_i =
      weights_vector.subvec(a, b) / variance_vector.subvec(a, b);
    arma::mat d_scaled_i = d_matrix_i;
    d_scaled_i.each_col() %= scale_i;
    ans += d_matrix_i.t() * d_scaled_i;
  }
  return ans / phi;
}
//==============================================================================


//============================ get_sc_criteria =================================
// [[Rcpp::export]]
Rcpp::List get_gee_criteria_sc_cw(const arma::vec& y_vector,
                                  const arma::vec& id_vector,
                                  const arma::vec& repeated_vector,
                                  const char* family,
                                  const arma::vec& mu_vector,
                                  const char* correlation_structure,
                                  const arma::vec& alpha_vector,
                                  const double& phi,
                                  const arma::vec& weights_vector) {
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  double sc_criterion = 0.0;
  double sum_log_det = 0.0;
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::vec s_vector_i = s_vector.subvec(a, b);
    const arma::mat v_matrix_i =
      get_v_matrix_cc(family,
                      mu_vector.subvec(a, b),
                      repeated_vector.subvec(a, b),
                      phi,
                      correlation_matrix,
                      weights_vector.subvec(a, b));
    sc_criterion +=
      arma::dot(s_vector_i, solve_chol_or_lu_vec(v_matrix_i, s_vector_i));
    double ld;
    double ld_sign;
    arma::log_det(ld, ld_sign, v_matrix_i);
    if (ld_sign <= 0.0) {
      Rcpp::stop(
        "get_gee_criteria_sc_cw: working covariance V_i is not positive "
        "definite for cluster at index %d -- check correlation parameters "
        "and model specification.",
        static_cast<int>(cl.start)
      );
    }
    sum_log_det += ld;
  }
  return Rcpp::List::create(
    Rcpp::Named("sc") = sc_criterion,
    Rcpp::Named("gp") = -0.5 * (sc_criterion + sum_log_det)
  );
}
//==============================================================================


//============================ sc criteria with odds ratios ====================
// [[Rcpp::export]]
Rcpp::List get_gee_criteria_sc_cw_or(const arma::vec& y_vector,
                                     const arma::vec& id_vector,
                                     const arma::vec& repeated_vector,
                                     const arma::vec& mu_vector,
                                     const arma::vec& alpha_vector,
                                     const arma::vec& weights_vector) {
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  double sc_criterion = 0.0;
  double sum_log_det = 0.0;
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::vec s_vector_i = s_vector.subvec(a, b);
    const arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector.subvec(a, b),
                                       repeated_max,
                                       alpha_vector);
    const arma::mat v_matrix_i =
      get_v_matrix_or(mu_vector.subvec(a, b),
                      odds_ratios_vector_i,
                      weights_vector.subvec(a, b));
    sc_criterion +=
      arma::dot(s_vector_i, solve_chol_or_lu_vec(v_matrix_i, s_vector_i));
    double ld;
    double ld_sign;
    arma::log_det(ld, ld_sign, v_matrix_i);
    if (ld_sign <= 0.0) {
      Rcpp::stop(
        "get_gee_criteria_sc_cw_or: working covariance V_i is not positive "
        "definite for cluster at index %d -- check odds-ratio parameters "
        "and model specification.",
        static_cast<int>(cl.start)
      );
    }
    sum_log_det += ld;
  }
  return Rcpp::List::create(
    Rcpp::Named("sc") = sc_criterion,
    Rcpp::Named("gp") = -0.5 * (sc_criterion + sum_log_det)
  );
}
//==============================================================================
