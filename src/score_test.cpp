#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link_functions.h"
#include "nuisance_quantities_cc.h"
#include "nuisance_quantities_or.h"
#include "clusterutils.h"
#include "utils.h"

//============================ estimating equations - cc =======================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec estimating_equations_gee_cc(const arma::vec& y_vector,
                                      const arma::mat& model_matrix,
                                      const arma::vec& id_vector,
                                      const arma::vec& repeated_vector,
                                      const arma::vec& weights_vector,
                                      const char* link,
                                      const char* family,
                                      const arma::vec& mu_vector,
                                      const arma::vec& eta_vector,
                                      const char* correlation_structure,
                                      const arma::vec& alpha_vector,
                                      const double& phi) {
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
  const arma::vec delta_vector = mueta(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::vec ans(params_no, arma::fill::zeros);
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
    const arma::mat v_matrix_i =
      get_v_matrix_cc(family,
                      mu_vector.subvec(a, b),
                      repeated_vector.subvec(a, b),
                      phi,
                      correlation_matrix,
                      weights_vector.subvec(a, b));
    ans += d_matrix_i.t() *
      solve_chol_or_lu_vec(v_matrix_i, s_vector.subvec(a, b));
  }
  return ans;
}
//==============================================================================


//============================ estimating equations - or =======================
// [[Rcpp::export]]
arma::vec estimating_equations_gee_or(const arma::vec& y_vector,
                                      const arma::mat& model_matrix,
                                      const arma::vec& id_vector,
                                      const arma::vec& repeated_vector,
                                      const arma::vec& weights_vector,
                                      const char* link,
                                      const arma::vec& mu_vector,
                                      const arma::vec& eta_vector,
                                      const arma::vec& alpha_vector) {
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  const arma::vec delta_vector = mueta(link, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::vec ans(params_no, arma::fill::zeros);
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
    const arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector.subvec(a, b),
                                       repeated_max,
                                       alpha_vector);
    const arma::mat v_matrix_i =
      get_v_matrix_or(mu_vector.subvec(a, b),
                      odds_ratios_vector_i,
                      weights_vector.subvec(a, b));
    ans += d_matrix_i.t() *
      solve_chol_or_lu_vec(v_matrix_i, s_vector.subvec(a, b));
  }
  return ans;
}
//==============================================================================
