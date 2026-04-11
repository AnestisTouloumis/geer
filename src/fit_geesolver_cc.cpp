#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "nuisance_quantities_cc.h"
#include "utils.h"
#include "link_functions.h"
#include "link_utils.h"
#include "family_utils.h"
#include "variance_functions.h"
#include "covariance_matrices.h"
#include "clusterutils.h"
#include "method_dispatch.h"
#include <algorithm>
#include <cmath>


//============================ update beta - gee ===============================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_beta_gee_cc(const arma::vec& y_vector,
                             const arma::mat& model_matrix,
                             const arma::vec& id_vector,
                             const arma::vec& repeated_vector,
                             const arma::vec& weights_vector,
                             const char* link,
                             const char* family,
                             const arma::vec& beta_vector,
                             const arma::vec& mu_vector,
                             const arma::vec& eta_vector,
                             const char* correlation_structure,
                             const arma::vec& alpha_vector,
                             const double& phi) {
  const LinkCode lc = parse_link(link);
  const FamilyCode fc = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat naive_matrix_inverse(params_no, params_no, arma::fill::zeros);
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
  const arma::vec delta_vector = mueta(lc, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::mat d_matrix_i;
  arma::mat v_matrix_i;
  arma::mat v_matrix_inverse_d_matrix_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::uword m = cl.end - cl.start;
    const arma::subview_col<double> delta_vector_i = delta_vector.subvec(a, b);
    const arma::subview_col<double> s_vector_i = s_vector.subvec(a, b);
    if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
      d_matrix_i.set_size(m, params_no);
    }
    d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector_i;
    v_matrix_i = get_v_matrix_cc(fc,
                                 mu_vector.subvec(a, b),
                                 repeated_vector.subvec(a, b),
                                 phi,
                                 correlation_matrix,
                                 weights_vector.subvec(a, b));
    if (v_matrix_inverse_d_matrix_i.n_rows != m ||
        v_matrix_inverse_d_matrix_i.n_cols != params_no) {
      v_matrix_inverse_d_matrix_i.set_size(m, params_no);
    }
    const bool ok_cluster = arma::solve(v_matrix_inverse_d_matrix_i,
                                        v_matrix_i,
                                        d_matrix_i,
                                        arma::solve_opts::likely_sympd);
    if (!ok_cluster) {
      Rcpp::stop("update_beta_gee_cc: failed to solve V_i^{-1} D_i.");
    }
    if (d_matrix_trans_v_matrix_inverse_i.n_rows != params_no ||
        d_matrix_trans_v_matrix_inverse_i.n_cols != m) {
      d_matrix_trans_v_matrix_inverse_i.set_size(params_no, m);
    }
    d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector_i;
  }
  symmetrize_if_close(naive_matrix_inverse, 1e-10);
  arma::vec update_step(params_no, arma::fill::zeros);
  const bool ok_step = arma::solve(update_step,
                                   naive_matrix_inverse,
                                   u_vector,
                                   arma::solve_opts::likely_sympd);
  if (!ok_step) {
    Rcpp::stop("update_beta_gee_cc: failed to solve for the beta update.");
  }
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - naive =============================
// [[Rcpp::export]]
arma::vec update_beta_naive_cc(const arma::vec& y_vector,
                               const arma::mat& model_matrix,
                               const arma::vec& id_vector,
                               const arma::vec& repeated_vector,
                               const arma::vec& weights_vector,
                               const char* link,
                               const char* family,
                               const arma::vec& beta_vector,
                               const arma::vec& mu_vector,
                               const arma::vec& eta_vector,
                               const char* correlation_structure,
                               const arma::vec& alpha_vector,
                               const double& phi) {
  const LinkCode lc = parse_link(link);
  const FamilyCode fc = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat lambda_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat naive_matrix_inverse(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = mueta(lc, eta_vector);
  const arma::vec delta_star_vector =
    mueta2(lc, eta_vector) / arma::square(delta_vector);
  const arma::vec variance_vector = variance(fc, mu_vector);
  const arma::vec alpha_star_vector =
    -0.5 * variancemu(fc, mu_vector) / variance_vector;
    const arma::vec alpha_star_plus_delta_star_vector =
    alpha_star_vector + delta_star_vector;
    const arma::vec s_vector = y_vector - mu_vector;
    const arma::mat correlation_matrix =
      get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
    const auto clusters = clusters_from_sorted_id(id_vector);
    arma::mat d_matrix_i;
    arma::mat v_matrix_i;
    arma::mat v_matrix_inverse_d_matrix_i;
    arma::mat d_matrix_trans_v_matrix_inverse_i;
    arma::mat v_matrix_inverse_alpha_plus_delta_star_diag_i;
    arma::mat v_matrix_inverse_delta_star_diag_i;
    arma::mat kron_d_matrix_d_matrix_i;
    arma::mat alpha_star_plus_delta_star_matrix_i;
    arma::mat delta_star_matrix_i;
    for (const auto& cl : clusters) {
      const arma::uword a = cl.start;
      const arma::uword b = cl.end - 1;
      const arma::uword m = cl.end - cl.start;
      const arma::subview_col<double> delta_vector_i = delta_vector.subvec(a, b);
      const arma::subview_col<double> s_vector_i = s_vector.subvec(a, b);
      const arma::subview_col<double> delta_star_vector_i =
        delta_star_vector.subvec(a, b);
      const arma::subview_col<double> alpha_star_plus_delta_star_vector_i =
        alpha_star_plus_delta_star_vector.subvec(a, b);
      if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
        d_matrix_i.set_size(m, params_no);
      }
      d_matrix_i = model_matrix.rows(a, b);
      d_matrix_i.each_col() %= delta_vector_i;
      v_matrix_i = get_v_matrix_cc(fc,
                                   mu_vector.subvec(a, b),
                                   repeated_vector.subvec(a, b),
                                   phi,
                                   correlation_matrix,
                                   weights_vector.subvec(a, b));

      if (v_matrix_inverse_d_matrix_i.n_rows != m ||
          v_matrix_inverse_d_matrix_i.n_cols != params_no) {
        v_matrix_inverse_d_matrix_i.set_size(m, params_no);
      }
      const bool ok_d = arma::solve(v_matrix_inverse_d_matrix_i,
                                    v_matrix_i,
                                    d_matrix_i,
                                    arma::solve_opts::likely_sympd);
      if (!ok_d) {
        Rcpp::stop("update_beta_naive_cc: failed to solve V_i^{-1} D_i.");
      }
      if (d_matrix_trans_v_matrix_inverse_i.n_rows != params_no ||
          d_matrix_trans_v_matrix_inverse_i.n_cols != m) {
        d_matrix_trans_v_matrix_inverse_i.set_size(params_no, m);
      }
      d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
      naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector_i;
      if (alpha_star_plus_delta_star_matrix_i.n_rows != m ||
          alpha_star_plus_delta_star_matrix_i.n_cols != m) {
        alpha_star_plus_delta_star_matrix_i.set_size(m, m);
      }
      alpha_star_plus_delta_star_matrix_i.zeros();
      alpha_star_plus_delta_star_matrix_i.diag() =
        alpha_star_plus_delta_star_vector_i;
      if (delta_star_matrix_i.n_rows != m || delta_star_matrix_i.n_cols != m) {
        delta_star_matrix_i.set_size(m, m);
      }
      delta_star_matrix_i.zeros();
      delta_star_matrix_i.diag() = delta_star_vector_i;
      if (v_matrix_inverse_alpha_plus_delta_star_diag_i.n_rows != m ||
          v_matrix_inverse_alpha_plus_delta_star_diag_i.n_cols != m) {
        v_matrix_inverse_alpha_plus_delta_star_diag_i.set_size(m, m);
      }
      const bool ok_alpha = arma::solve(v_matrix_inverse_alpha_plus_delta_star_diag_i,
                                        v_matrix_i,
                                        alpha_star_plus_delta_star_matrix_i,
                                        arma::solve_opts::likely_sympd);
      if (!ok_alpha) {
        Rcpp::stop("update_beta_naive_cc: failed to solve V_i^{-1} diag(alpha* + delta*).");
      }
      if (v_matrix_inverse_delta_star_diag_i.n_rows != m ||
          v_matrix_inverse_delta_star_diag_i.n_cols != m) {
        v_matrix_inverse_delta_star_diag_i.set_size(m, m);
      }
      const bool ok_delta = arma::solve(v_matrix_inverse_delta_star_diag_i,
                                        v_matrix_i,
                                        delta_star_matrix_i,
                                        arma::solve_opts::likely_sympd);
      if (!ok_delta) {
        Rcpp::stop("update_beta_naive_cc: failed to solve V_i^{-1} diag(delta*).");
      }
      kron_self_matrix_into(kron_d_matrix_d_matrix_i, d_matrix_i);
      lambda_matrix +=
        kron_d_matrix_d_matrix_i.t() *
        (kappa_right(v_matrix_inverse_alpha_plus_delta_star_diag_i.t()) -
        kronecker_identity_right_kappa(v_matrix_inverse_alpha_plus_delta_star_diag_i) -
        kronecker_left_identity_kappa(v_matrix_inverse_delta_star_diag_i)) *
        d_matrix_i;
    }
    symmetrize_if_close(naive_matrix_inverse, 1e-10);
    const arma::vec lambda_vector =
      lambda_from_blocks_chol_or_lu(naive_matrix_inverse, lambda_matrix);
    arma::vec update_step(params_no, arma::fill::zeros);
    const bool ok_step = arma::solve(update_step,
                                     naive_matrix_inverse,
                                     u_vector + lambda_vector,
                                     arma::solve_opts::likely_sympd);
    if (!ok_step) {
      Rcpp::stop("update_beta_naive_cc: failed to solve for the beta update.");
    }
    return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - robust  ===========================
// [[Rcpp::export]]
arma::vec update_beta_robust_cc(const arma::vec& y_vector,
                                const arma::mat& model_matrix,
                                const arma::vec& id_vector,
                                const arma::vec& repeated_vector,
                                const arma::vec& weights_vector,
                                const char* link,
                                const char* family,
                                const arma::vec& beta_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& eta_vector,
                                const char* correlation_structure,
                                const arma::vec& alpha_vector,
                                const double& phi) {
  const LinkCode lc = parse_link(link);
  const FamilyCode fc = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat partial_derivatives_matrix(params_no * params_no,
                                       params_no,
                                       arma::fill::zeros);
  arma::mat second_derivatives_matrix(params_no * params_no,
                                      params_no,
                                      arma::fill::zeros);
  arma::mat naive_matrix_inverse(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = mueta(lc, eta_vector);
  const arma::vec delta_star_vector =
    mueta2(lc, eta_vector) / arma::square(delta_vector);
  const arma::vec variance_vector = variance(fc, mu_vector);
  const arma::vec alpha_star_vector =
    -0.5 * variancemu(fc, mu_vector) / variance_vector;
    const arma::vec alpha_star_plus_delta_star_vector =
    alpha_star_vector + delta_star_vector;
    const arma::vec s_vector = y_vector - mu_vector;
    const arma::mat correlation_matrix =
      get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
    const auto clusters = clusters_from_sorted_id(id_vector);
    arma::mat d_matrix_i;
    arma::mat v_matrix_i;
    arma::mat v_matrix_inverse_d_matrix_i;
    arma::mat d_matrix_trans_v_matrix_inverse_i;
    arma::vec u_vector_i(params_no, arma::fill::zeros);
    arma::mat alpha_star_plus_delta_star_matrix_i;
    arma::mat alpha_star_matrix_i;
    arma::mat v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i;
    arma::mat v_matrix_inverse_alpha_star_matrix_i;
    arma::mat kron_d_matrix_d_matrix_i;
    arma::mat kron_d_matrix_trans_d_matrix_trans_i;
    for (const auto& cl : clusters) {
      const arma::uword a = cl.start;
      const arma::uword b = cl.end - 1;
      const arma::uword m = cl.end - cl.start;
      const arma::subview_col<double> delta_vector_i = delta_vector.subvec(a, b);
      const arma::subview_col<double> s_vector_i = s_vector.subvec(a, b);
      const arma::subview_col<double> alpha_star_vector_i =
        alpha_star_vector.subvec(a, b);
      const arma::subview_col<double> delta_star_vector_i =
        delta_star_vector.subvec(a, b);
      const arma::subview_col<double> alpha_star_plus_delta_star_vector_i =
        alpha_star_plus_delta_star_vector.subvec(a, b);
      if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
        d_matrix_i.set_size(m, params_no);
      }
      d_matrix_i = model_matrix.rows(a, b);
      d_matrix_i.each_col() %= delta_vector_i;
      v_matrix_i = get_v_matrix_cc(fc,
                                   mu_vector.subvec(a, b),
                                   repeated_vector.subvec(a, b),
                                   phi,
                                   correlation_matrix,
                                   weights_vector.subvec(a, b));
      if (v_matrix_inverse_d_matrix_i.n_rows != m ||
          v_matrix_inverse_d_matrix_i.n_cols != params_no) {
        v_matrix_inverse_d_matrix_i.set_size(m, params_no);
      }
      const bool ok_d = arma::solve(v_matrix_inverse_d_matrix_i,
                                    v_matrix_i,
                                    d_matrix_i,
                                    arma::solve_opts::likely_sympd);
      if (!ok_d) {
        Rcpp::stop("update_beta_robust_cc: failed to solve V_i^{-1} D_i.");
      }
      if (d_matrix_trans_v_matrix_inverse_i.n_rows != params_no ||
          d_matrix_trans_v_matrix_inverse_i.n_cols != m) {
        d_matrix_trans_v_matrix_inverse_i.set_size(params_no, m);
      }
      d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
      naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector_i;
      u_vector += u_vector_i;
      meat_matrix += u_vector_i * u_vector_i.t();
      if (alpha_star_plus_delta_star_matrix_i.n_rows != m ||
          alpha_star_plus_delta_star_matrix_i.n_cols != m) {
        alpha_star_plus_delta_star_matrix_i.set_size(m, m);
      }
      alpha_star_plus_delta_star_matrix_i.zeros();
      alpha_star_plus_delta_star_matrix_i.diag() =
        alpha_star_plus_delta_star_vector_i;
      if (alpha_star_matrix_i.n_rows != m || alpha_star_matrix_i.n_cols != m) {
        alpha_star_matrix_i.set_size(m, m);
      }
      alpha_star_matrix_i.zeros();
      alpha_star_matrix_i.diag() = alpha_star_vector_i;
      if (v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i.n_rows != m ||
          v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i.n_cols != m) {
        v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i.set_size(m, m);
      }
      const bool ok_alpha_plus_delta =
        arma::solve(v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i,
                    v_matrix_i,
                    alpha_star_plus_delta_star_matrix_i,
                    arma::solve_opts::likely_sympd);
      if (!ok_alpha_plus_delta) {
        Rcpp::stop(
          "update_beta_robust_cc: failed to solve V_i^{-1} diag(alpha* + delta*)."
        );
      }
      if (v_matrix_inverse_alpha_star_matrix_i.n_rows != m ||
          v_matrix_inverse_alpha_star_matrix_i.n_cols != m) {
        v_matrix_inverse_alpha_star_matrix_i.set_size(m, m);
      }
      const bool ok_alpha =
        arma::solve(v_matrix_inverse_alpha_star_matrix_i,
                    v_matrix_i,
                    alpha_star_matrix_i,
                    arma::solve_opts::likely_sympd);
      if (!ok_alpha) {
        Rcpp::stop("update_beta_robust_cc: failed to solve V_i^{-1} diag(alpha*).");
      }
      const arma::mat
      kappa_matrix_delta_star_matrix_plus_alpha_star_matrix_v_matrix_inverse_i =
        kappa_right(v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i.t());
      kron_self_matrix_into(kron_d_matrix_d_matrix_i, d_matrix_i);
      if (kron_d_matrix_trans_d_matrix_trans_i.n_rows !=
          kron_d_matrix_d_matrix_i.n_cols ||
          kron_d_matrix_trans_d_matrix_trans_i.n_cols !=
          kron_d_matrix_d_matrix_i.n_rows) {
        kron_d_matrix_trans_d_matrix_trans_i.set_size(
          kron_d_matrix_d_matrix_i.n_cols,
          kron_d_matrix_d_matrix_i.n_rows
        );
      }
      kron_d_matrix_trans_d_matrix_trans_i = kron_d_matrix_d_matrix_i.t();
      second_derivatives_matrix -=
        kron_d_matrix_trans_d_matrix_trans_i *
        (kappa_matrix_delta_star_matrix_plus_alpha_star_matrix_v_matrix_inverse_i +
        kronecker_left_identity_kappa(
          v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i +
            v_matrix_inverse_alpha_star_matrix_i
        ) +
          kronecker_identity_right_kappa(
            v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i
          )) *
            d_matrix_i;
      partial_derivatives_matrix +=
        kron_d_matrix_trans_d_matrix_trans_i *
        (kappa_matrix_delta_star_matrix_plus_alpha_star_matrix_v_matrix_inverse_i +
        kronecker_left_identity_kappa(v_matrix_inverse_alpha_star_matrix_i)) *
        s_vector_i * u_vector_i.t();
    }
    symmetrize_if_close(naive_matrix_inverse, 1e-10);
    arma::mat naive_matrix_meat_matrix(params_no, params_no, arma::fill::zeros);
    const bool ok_naive_meat = arma::solve(naive_matrix_meat_matrix,
                                           naive_matrix_inverse,
                                           meat_matrix,
                                           arma::solve_opts::likely_sympd);
    if (!ok_naive_meat) {
      Rcpp::stop("update_beta_robust_cc: failed to solve bread^{-1} meat.");
    }
    arma::mat robust_matrix(params_no, params_no, arma::fill::zeros);
    const bool ok_robust = arma::solve(robust_matrix,
                                       naive_matrix_inverse,
                                       naive_matrix_meat_matrix.t(),
                                       arma::solve_opts::likely_sympd);
    if (!ok_robust) {
      Rcpp::stop("update_beta_robust_cc: failed to solve for robust covariance.");
    }
    arma::vec lambda_vector(params_no, arma::fill::zeros);
    for (arma::uword r = 0; r < params_no; ++r) {
      const arma::mat first_block =
        partial_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
      const arma::mat second_block =
        second_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
      arma::mat naive_matrix_first_block(params_no, params_no, arma::fill::zeros);
      const bool ok_first_block = arma::solve(naive_matrix_first_block,
                                              naive_matrix_inverse,
                                              first_block,
                                              arma::solve_opts::likely_sympd);
      if (!ok_first_block) {
        Rcpp::stop("update_beta_robust_cc: failed to solve bread^{-1} first block.");
      }
      lambda_vector[r] =
        -(arma::trace(naive_matrix_first_block) +
        0.5 * arma::trace(robust_matrix * second_block));
    }
    arma::vec update_step(params_no, arma::fill::zeros);
    const bool ok_step = arma::solve(update_step,
                                     naive_matrix_inverse,
                                     u_vector + lambda_vector,
                                     arma::solve_opts::likely_sympd);
    if (!ok_step) {
      Rcpp::stop("update_beta_robust_cc: failed to solve for the beta update.");
    }
    return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - empirical =========================
// [[Rcpp::export]]
arma::vec update_beta_empirical_cc(const arma::vec& y_vector,
                                   const arma::mat& model_matrix,
                                   const arma::vec& id_vector,
                                   const arma::vec& repeated_vector,
                                   const arma::vec& weights_vector,
                                   const char* link,
                                   const char* family,
                                   const arma::vec& beta_vector,
                                   const arma::vec& mu_vector,
                                   const arma::vec& eta_vector,
                                   const char* correlation_structure,
                                   const arma::vec& alpha_vector,
                                   const double& phi) {
  const LinkCode lc = parse_link(link);
  const FamilyCode fc = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat partial_derivatives_matrix(params_no * params_no,
                                       params_no,
                                       arma::fill::zeros);
  arma::mat second_derivatives_matrix(params_no * params_no,
                                      params_no,
                                      arma::fill::zeros);
  arma::mat observed_fisher_info_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec s_vector = y_vector - mu_vector;
  const arma::vec delta_vector = mueta(lc, eta_vector);
  const arma::vec mueta2_vector = mueta2(lc, eta_vector);
  const arma::vec delta_star_vector =
    mueta2_vector / arma::square(delta_vector);
  const arma::vec delta_tilde_star_vector =
    (delta_vector % mueta3(lc, eta_vector) -
    2.0 * arma::square(mueta2_vector)) /
      arma::square(arma::square(delta_vector));
  const arma::vec variance_vector = variance(fc, mu_vector);
  const arma::vec variancemu_vector = variancemu(fc, mu_vector);
  const arma::vec alpha_star_vector =
    -0.5 * variancemu_vector / variance_vector;
    const arma::vec alpha_tilde_star_vector =
    0.5 * (arma::square(variancemu_vector) / variance_vector -
    variancemu2(fc, mu_vector)) / variance_vector;
    const arma::vec alpha_star_plus_delta_star_vector =
      alpha_star_vector + delta_star_vector;
    const arma::mat correlation_matrix =
      get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
    const auto clusters = clusters_from_sorted_id(id_vector);
    arma::mat d_matrix_i;
    arma::mat v_matrix_i;
    arma::mat v_matrix_inverse_d_matrix_i;
    arma::mat v_matrix_inverse_i;
    arma::mat epsilon_matrix_i;
    arma::mat observed_fisher_info_matrix_i;
    arma::mat alpha_star_plus_delta_star_matrix_i;
    arma::mat temp_diagonal_matrix;
    arma::mat kron_d_matrix_d_matrix_i;
    arma::mat kron_d_matrix_trans_d_matrix_trans_i;
    arma::mat hessian_term2;
    arma::mat hessian_term3;
    arma::mat diagonal_correction;
    arma::vec u_vector_i(params_no, arma::fill::zeros);
    arma::vec weighted_residuals_i;
    for (const auto& cl : clusters) {
      const arma::uword a = cl.start;
      const arma::uword b = cl.end - 1;
      const arma::uword m = cl.end - cl.start;
      const arma::vec s_vector_i = s_vector.subvec(a, b);
      const arma::vec alpha_star_vector_i = alpha_star_vector.subvec(a, b);
      const arma::vec alpha_tilde_star_vector_i =
        alpha_tilde_star_vector.subvec(a, b);
      const arma::vec alpha_star_plus_delta_star_vector_i =
        alpha_star_plus_delta_star_vector.subvec(a, b);
      if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
        d_matrix_i.set_size(m, params_no);
      }
      d_matrix_i = model_matrix.rows(a, b);
      d_matrix_i.each_col() %= delta_vector.subvec(a, b);

      v_matrix_i = get_v_matrix_cc(fc,
                                   mu_vector.subvec(a, b),
                                   repeated_vector.subvec(a, b),
                                   phi,
                                   correlation_matrix,
                                   weights_vector.subvec(a, b));
      if (v_matrix_inverse_d_matrix_i.n_rows != m ||
          v_matrix_inverse_d_matrix_i.n_cols != params_no) {
        v_matrix_inverse_d_matrix_i.set_size(m, params_no);
      }
      const bool ok_d = arma::solve(v_matrix_inverse_d_matrix_i,
                                    v_matrix_i,
                                    d_matrix_i,
                                    arma::solve_opts::likely_sympd);
      if (!ok_d) {
        Rcpp::stop("update_beta_empirical_cc: failed to solve V_i^{-1} D_i.");
      }
      u_vector_i = v_matrix_inverse_d_matrix_i.t() * s_vector_i;
      u_vector += u_vector_i;
      meat_matrix += u_vector_i * u_vector_i.t();
      v_matrix_inverse_i.set_size(m, m);
      const bool ok_vinv = arma::solve(v_matrix_inverse_i,
                                       v_matrix_i,
                                       arma::eye(m, m),
                                       arma::solve_opts::likely_sympd);
      if (!ok_vinv) {
        Rcpp::stop("update_beta_empirical_cc: failed to solve V_i^{-1}.");
      }
      weighted_residuals_i = v_matrix_inverse_i * s_vector_i;
      epsilon_matrix_i.set_size(m, m);
      epsilon_matrix_i.zeros();
      epsilon_matrix_i.diag() =
        weighted_residuals_i % alpha_star_plus_delta_star_vector_i;
      temp_diagonal_matrix.set_size(m, m);
      temp_diagonal_matrix.zeros();
      temp_diagonal_matrix.diag() = s_vector_i % alpha_star_vector_i - 1.0;
      epsilon_matrix_i += v_matrix_inverse_i * temp_diagonal_matrix;
      observed_fisher_info_matrix_i = d_matrix_i.t() * epsilon_matrix_i * d_matrix_i;
      observed_fisher_info_matrix -= observed_fisher_info_matrix_i;
      partial_derivatives_matrix +=
        arma::vectorise(observed_fisher_info_matrix_i.t()) * u_vector_i.t();
      alpha_star_plus_delta_star_matrix_i.set_size(m, m);
      alpha_star_plus_delta_star_matrix_i.zeros();
      alpha_star_plus_delta_star_matrix_i.diag() =
        alpha_star_plus_delta_star_vector_i;
      kron_self_matrix_into(kron_d_matrix_d_matrix_i, d_matrix_i);
      kron_d_matrix_trans_d_matrix_trans_i = kron_d_matrix_d_matrix_i.t();
      const arma::mat hessian_term1 =
        alpha_star_plus_delta_star_matrix_i * epsilon_matrix_i;
      const arma::vec hessian_diag2 =
        alpha_star_plus_delta_star_vector_i %
        (alpha_star_plus_delta_star_vector_i + alpha_star_vector_i) %
        weighted_residuals_i;
      hessian_term2.set_size(m, m);
      hessian_term2.zeros();
      hessian_term2.diag() = -hessian_diag2;
      const arma::vec hessian_diag3 =
        (alpha_tilde_star_vector_i + delta_tilde_star_vector.subvec(a, b)) %
        weighted_residuals_i;
      hessian_term3.set_size(m, m);
      hessian_term3.zeros();
      hessian_term3.diag() = hessian_diag3;
      const arma::mat kappa_correction =
        kappa_right(hessian_term1 + hessian_term2 + hessian_term3);
      const arma::mat right_kronecker_correction =
        kronecker_identity_right_kappa(
          epsilon_matrix_i.t() * alpha_star_plus_delta_star_matrix_i
        );
      diagonal_correction.set_size(m, m);
      diagonal_correction.zeros();
      diagonal_correction.diag() =
        s_vector_i % alpha_tilde_star_vector_i - alpha_star_vector_i;
      const arma::mat left_kronecker_correction =
        kronecker_left_identity_kappa(
          epsilon_matrix_i * alpha_star_plus_delta_star_matrix_i +
            v_matrix_inverse_i * diagonal_correction
        );
      second_derivatives_matrix +=
        kron_d_matrix_trans_d_matrix_trans_i *
        (kappa_correction +
        right_kronecker_correction +
        left_kronecker_correction) *
        d_matrix_i;
    }
    arma::mat robust_matrix(params_no, params_no, arma::fill::zeros);
    arma::vec lambda_vector(params_no, arma::fill::zeros);
    arma::mat lu_lower, lu_upper, permutation_matrix;
    const bool lu_success =
      arma::lu(lu_lower, lu_upper, permutation_matrix, observed_fisher_info_matrix);
    if (!lu_success) {
      Rcpp::stop("update_beta_empirical_cc: LU factorization failed.");
    }
    auto solve_observed_fisher_matrix = [&](const arma::mat& rhs_matrix) -> arma::mat {
      const arma::mat permuted_rhs = permutation_matrix * rhs_matrix;
      arma::mat lu_forward;
      const bool ok_forward =
        arma::solve(lu_forward,
                    arma::trimatl(lu_lower),
                    permuted_rhs,
                    arma::solve_opts::allow_ugly);
      if (!ok_forward) {
        Rcpp::stop("update_beta_empirical_cc: forward LU solve failed.");
      }
      arma::mat ans;
      const bool ok_back =
        arma::solve(ans,
                    arma::trimatu(lu_upper),
                    lu_forward,
                    arma::solve_opts::allow_ugly);
      if (!ok_back) {
        Rcpp::stop("update_beta_empirical_cc: backward LU solve failed.");
      }
      return ans;
    };
    auto solve_observed_fisher_vector = [&](const arma::vec& rhs_vector) -> arma::vec {
      const arma::vec permuted_rhs = permutation_matrix * rhs_vector;
      arma::vec lu_forward;
      const bool ok_forward =
        arma::solve(lu_forward,
                    arma::trimatl(lu_lower),
                    permuted_rhs,
                    arma::solve_opts::allow_ugly);
      if (!ok_forward) {
        Rcpp::stop("update_beta_empirical_cc: forward LU solve failed.");
      }
      arma::vec ans;
      const bool ok_back =
        arma::solve(ans,
                    arma::trimatu(lu_upper),
                    lu_forward,
                    arma::solve_opts::allow_ugly);
      if (!ok_back) {
        Rcpp::stop("update_beta_empirical_cc: backward LU solve failed.");
      }
      return ans;
    };
    const arma::mat fisher_inv_meat = solve_observed_fisher_matrix(meat_matrix);
    robust_matrix = solve_observed_fisher_matrix(fisher_inv_meat.t());
    for (arma::uword param_idx = 0; param_idx < params_no; ++param_idx) {
      const arma::uword block_start = param_idx * params_no;
      const arma::uword block_end = (param_idx + 1) * params_no - 1;
      const arma::mat partial_deriv_block =
        partial_derivatives_matrix.rows(block_start, block_end);
      const arma::mat second_deriv_block =
        second_derivatives_matrix.rows(block_start, block_end);
      lambda_vector[param_idx] =
        -(arma::trace(solve_observed_fisher_matrix(partial_deriv_block)) +
        0.5 * arma::trace(robust_matrix * second_deriv_block));
    }
    const arma::vec update_step =
      solve_observed_fisher_vector(u_vector + lambda_vector);
    return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - jeffreys ==========================
// [[Rcpp::export]]
arma::vec update_beta_jeffreys_cc(const arma::vec& y_vector,
                                  const arma::mat& model_matrix,
                                  const arma::vec& id_vector,
                                  const arma::vec& repeated_vector,
                                  const arma::vec& weights_vector,
                                  const char* link,
                                  const char* family,
                                  const arma::vec& beta_vector,
                                  const arma::vec& mu_vector,
                                  const arma::vec& eta_vector,
                                  const char* correlation_structure,
                                  const arma::vec& alpha_vector,
                                  const double& phi,
                                  const double& jeffreys_power) {
  const LinkCode lc = parse_link(link);
  const FamilyCode fc = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat naive_matrix_inverse(params_no, params_no, arma::fill::zeros);
  arma::mat lambda_matrix(params_no * params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = mueta(lc, eta_vector);
  const arma::vec delta_star_vector =
    mueta2(lc, eta_vector) / arma::square(delta_vector);
  const arma::vec variance_vector = variance(fc, mu_vector);
  const arma::vec alpha_star_vector =
    -0.5 * variancemu(fc, mu_vector) / variance_vector;
    const arma::vec alpha_star_plus_delta_star_vector =
    alpha_star_vector + delta_star_vector;
    const arma::vec s_vector = y_vector - mu_vector;
    const arma::mat correlation_matrix =
      get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
    const auto clusters = clusters_from_sorted_id(id_vector);
    arma::mat d_matrix_i;
    arma::mat v_matrix_i;
    arma::mat v_matrix_inverse_d_matrix_i;
    arma::mat d_matrix_trans_v_matrix_inverse_i;
    arma::mat alpha_star_plus_delta_star_matrix_i;
    arma::mat v_matrix_inverse_alpha_star_plus_delta_star_matrix_i;
    arma::mat kron_d_matrix_d_matrix_i;
    arma::mat kron_d_matrix_trans_d_matrix_trans_i;
    for (const auto& cl : clusters) {
      const arma::uword a = cl.start;
      const arma::uword b = cl.end - 1;
      const arma::uword m = cl.end - cl.start;
      const arma::subview_col<double> delta_vector_i = delta_vector.subvec(a, b);
      const arma::subview_col<double> s_vector_i = s_vector.subvec(a, b);
      const arma::subview_col<double> alpha_star_plus_delta_star_vector_i =
        alpha_star_plus_delta_star_vector.subvec(a, b);

      if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
        d_matrix_i.set_size(m, params_no);
      }
      d_matrix_i = model_matrix.rows(a, b);
      d_matrix_i.each_col() %= delta_vector_i;
      v_matrix_i = get_v_matrix_cc(fc,
                                   mu_vector.subvec(a, b),
                                   repeated_vector.subvec(a, b),
                                   phi,
                                   correlation_matrix,
                                   weights_vector.subvec(a, b));
      if (v_matrix_inverse_d_matrix_i.n_rows != m ||
          v_matrix_inverse_d_matrix_i.n_cols != params_no) {
        v_matrix_inverse_d_matrix_i.set_size(m, params_no);
      }
      const bool ok_d = arma::solve(v_matrix_inverse_d_matrix_i,
                                    v_matrix_i,
                                    d_matrix_i,
                                    arma::solve_opts::likely_sympd);
      if (!ok_d) {
        Rcpp::stop("update_beta_jeffreys_cc: failed to solve V_i^{-1} D_i.");
      }
      if (d_matrix_trans_v_matrix_inverse_i.n_rows != params_no ||
          d_matrix_trans_v_matrix_inverse_i.n_cols != m) {
        d_matrix_trans_v_matrix_inverse_i.set_size(params_no, m);
      }
      d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
      naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector_i;
      if (alpha_star_plus_delta_star_matrix_i.n_rows != m ||
          alpha_star_plus_delta_star_matrix_i.n_cols != m) {
        alpha_star_plus_delta_star_matrix_i.set_size(m, m);
      }
      alpha_star_plus_delta_star_matrix_i.zeros();
      alpha_star_plus_delta_star_matrix_i.diag() =
        alpha_star_plus_delta_star_vector_i;
      if (v_matrix_inverse_alpha_star_plus_delta_star_matrix_i.n_rows != m ||
          v_matrix_inverse_alpha_star_plus_delta_star_matrix_i.n_cols != m) {
        v_matrix_inverse_alpha_star_plus_delta_star_matrix_i.set_size(m, m);
      }
      const bool ok_alpha = arma::solve(v_matrix_inverse_alpha_star_plus_delta_star_matrix_i,
                                        v_matrix_i,
                                        alpha_star_plus_delta_star_matrix_i,
                                        arma::solve_opts::likely_sympd);
      if (!ok_alpha) {
        Rcpp::stop(
          "update_beta_jeffreys_cc: failed to solve V_i^{-1} diag(alpha* + delta*)."
        );
      }
      kron_self_matrix_into(kron_d_matrix_d_matrix_i, d_matrix_i);
      if (kron_d_matrix_trans_d_matrix_trans_i.n_rows !=
          kron_d_matrix_d_matrix_i.n_cols ||
          kron_d_matrix_trans_d_matrix_trans_i.n_cols !=
          kron_d_matrix_d_matrix_i.n_rows) {
        kron_d_matrix_trans_d_matrix_trans_i.set_size(
          kron_d_matrix_d_matrix_i.n_cols,
          kron_d_matrix_d_matrix_i.n_rows
        );
      }
      kron_d_matrix_trans_d_matrix_trans_i = kron_d_matrix_d_matrix_i.t();
      lambda_matrix +=
        kron_d_matrix_trans_d_matrix_trans_i *
        (kronecker_left_identity_kappa(
            v_matrix_inverse_alpha_star_plus_delta_star_matrix_i
        ) +
          kronecker_identity_right_kappa(
            v_matrix_inverse_alpha_star_plus_delta_star_matrix_i
          )) *
            d_matrix_i;
    }
    symmetrize_if_close(naive_matrix_inverse, 1e-10);
    arma::mat naive_matrix(params_no, params_no, arma::fill::zeros);
    const bool ok_naive = arma::solve(naive_matrix,
                                      naive_matrix_inverse,
                                      arma::eye(params_no, params_no),
                                      arma::solve_opts::likely_sympd);
    if (!ok_naive) {
      Rcpp::stop("update_beta_jeffreys_cc: failed to invert naive_matrix_inverse.");
    }
    const arma::vec naive_matrix_vectorised = arma::vectorise(naive_matrix);
    const arma::vec lambda_vector =
      jeffreys_power * (lambda_matrix.t() * naive_matrix_vectorised);
    arma::vec update_step(params_no, arma::fill::zeros);
    const bool ok_step = arma::solve(update_step,
                                     naive_matrix_inverse,
                                     u_vector + lambda_vector,
                                     arma::solve_opts::likely_sympd);
    if (!ok_step) {
      Rcpp::stop("update_beta_jeffreys_cc: failed to solve for the beta update.");
    }
    return beta_vector + update_step;
}
//==============================================================================


//============================ update beta =====================================
// [[Rcpp::export]]
arma::vec update_beta_cc(const arma::vec& y_vector,
                         const arma::mat& model_matrix,
                         const arma::vec& id_vector,
                         const arma::vec& repeated_vector,
                         const arma::vec& weights_vector,
                         const char* link,
                         const char* family,
                         const arma::vec& beta_vector,
                         const arma::vec& mu_vector,
                         const arma::vec& eta_vector,
                         const char* correlation_structure,
                         const arma::vec& alpha_vector,
                         const double& phi,
                         const double& jeffreys_power,
                         const char* method) {
  switch (method_code(method)) {
  case M_GEE:
    return update_beta_gee_cc(y_vector,
                              model_matrix,
                              id_vector,
                              repeated_vector,
                              weights_vector,
                              link,
                              family,
                              beta_vector,
                              mu_vector,
                              eta_vector,
                              correlation_structure,
                              alpha_vector,
                              phi);
  case M_BR_NAIVE:
    return update_beta_naive_cc(y_vector,
                                model_matrix,
                                id_vector,
                                repeated_vector,
                                weights_vector,
                                link,
                                family,
                                beta_vector,
                                mu_vector,
                                eta_vector,
                                correlation_structure,
                                alpha_vector,
                                phi);
  case M_BR_ROBUST:
    return update_beta_robust_cc(y_vector,
                                 model_matrix,
                                 id_vector,
                                 repeated_vector,
                                 weights_vector,
                                 link,
                                 family,
                                 beta_vector,
                                 mu_vector,
                                 eta_vector,
                                 correlation_structure,
                                 alpha_vector,
                                 phi);
  case M_BR_EMPIRICAL:
    return update_beta_empirical_cc(y_vector,
                                    model_matrix,
                                    id_vector,
                                    repeated_vector,
                                    weights_vector,
                                    link,
                                    family,
                                    beta_vector,
                                    mu_vector,
                                    eta_vector,
                                    correlation_structure,
                                    alpha_vector,
                                    phi);
  case M_PGEE_JEFFREYS:
    return update_beta_jeffreys_cc(y_vector,
                                   model_matrix,
                                   id_vector,
                                   repeated_vector,
                                   weights_vector,
                                   link,
                                   family,
                                   beta_vector,
                                   mu_vector,
                                   eta_vector,
                                   correlation_structure,
                                   alpha_vector,
                                   phi,
                                   jeffreys_power);
  }
  Rcpp::stop("update_beta_cc: unsupported method.");
}
//==============================================================================


//=========================== fitting function =================================
// [[Rcpp::export]]
Rcpp::List fit_geesolver_cc(const arma::vec& y_vector,
                            const arma::mat& model_matrix,
                            const arma::vec& id_vector,
                            const arma::vec& repeated_vector,
                            const arma::vec& weights_vector,
                            const char* link,
                            const char* family,
                            arma::vec beta_vector,
                            const arma::vec& offset,
                            const int& maxiter,
                            const double& tolerance,
                            const int& step_maxiter,
                            const int& step_multiplier,
                            const double& jeffreys_power,
                            const char* method,
                            int use_params,
                            arma::vec alpha_vector,
                            const int& alpha_fixed,
                            const char* correlation_structure,
                            const int& mdependence,
                            double phi,
                            const int& phi_fixed) {
  const LinkCode lc = parse_link(link);
  const FamilyCode fc = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  use_params *= static_cast<int>(params_no);
  arma::vec beta_vector_new = beta_vector;
  arma::mat beta_hat_matrix = beta_vector;
  arma::vec stepsize_vector(params_no, arma::fill::zeros);
  arma::vec criterion_vector(maxiter, arma::fill::zeros);
  arma::vec beta_vector_inner(params_no, arma::fill::zeros);
  arma::vec beta_vector_new_inner(params_no, arma::fill::zeros);
  arma::vec stepsize_vector_inner(params_no, arma::fill::zeros);
  arma::vec eta_vector = model_matrix * beta_vector + offset;
  if (!valideta(link, arma2vec(eta_vector))) {
    Rcpp::stop(
      "invalid initial linear predictor: please another set of initial values for beta!!"
    );
  }
  arma::vec mu_vector = linkinv(lc, eta_vector);
  if (!validmu(family, arma2vec(mu_vector))) {
    Rcpp::stop(
      "invalid initial fitted values: please another set of initial values for beta!!"
    );
  }
  arma::vec pearson_residuals_vector =
    get_pearson_residuals(fc,
                          y_vector,
                          mu_vector,
                          weights_vector);

  if (phi_fixed == 0) {
    phi = get_phi_hat(pearson_residuals_vector, use_params);
  }
  if (alpha_fixed == 0) {
    alpha_vector = get_alpha_hat(correlation_structure,
                                 pearson_residuals_vector,
                                 id_vector,
                                 repeated_vector,
                                 phi,
                                 use_params,
                                 mdependence);
  }
  for (int i = 1; i < maxiter + 1; ++i) {
    stepsize_vector =
      update_beta_cc(y_vector,
                     model_matrix,
                     id_vector,
                     repeated_vector,
                     weights_vector,
                     link,
                     family,
                     beta_vector,
                     mu_vector,
                     eta_vector,
                     correlation_structure,
                     alpha_vector,
                     phi,
                     jeffreys_power,
                     method) - beta_vector;
    // Step-length control: monotone step-size rule (Ortega & Rheinboldt).
    // Accept the halved step if the resulting Newton step norm is strictly
    // smaller than the norm at the original iterate.  GEE has no log-
    // likelihood, so a merit-function (Armijo) line search is not available;
    // this norm-reduction criterion is the standard substitute used in geepack
    // and glmtoolbox.  criterion_inner is intentionally held fixed at the
    // original step norm throughout the inner loop — every candidate is
    // compared against that baseline, not against the previous candidate.
    double criterion_inner = arma::norm(stepsize_vector, "inf");
    beta_vector_inner = beta_vector;
    beta_vector_new = beta_vector;
    stepsize_vector_inner = stepsize_vector;
    for (int j = 1; j < step_maxiter + 1; ++j) {
      beta_vector_new_inner =
        beta_vector_inner +
        step_multiplier * std::pow(0.5, j - 1) * stepsize_vector_inner;
      eta_vector = model_matrix * beta_vector_new_inner + offset;
      if (!valideta(link, arma2vec(eta_vector))) {
        Rcpp::stop(
          "invalid linear predictor - please another set of initial values for beta!!"
        );
      }
      mu_vector = linkinv(lc, eta_vector);
      if (!validmu(family, arma2vec(mu_vector))) {
        Rcpp::stop(
          "invalid fitted values - please another set of initial values for beta!!"
        );
      }
      pearson_residuals_vector =
        get_pearson_residuals(fc,
                              y_vector,
                              mu_vector,
                              weights_vector);
      if (phi_fixed == 0) {
        phi = get_phi_hat(pearson_residuals_vector, use_params);
      }
      if (alpha_fixed == 0) {
        alpha_vector = get_alpha_hat(correlation_structure,
                                     pearson_residuals_vector,
                                     id_vector,
                                     repeated_vector,
                                     phi,
                                     use_params,
                                     mdependence);
      }
      stepsize_vector_inner =
        update_beta_cc(y_vector,
                       model_matrix,
                       id_vector,
                       repeated_vector,
                       weights_vector,
                       link,
                       family,
                       beta_vector_new_inner,
                       mu_vector,
                       eta_vector,
                       correlation_structure,
                       alpha_vector,
                       phi,
                       jeffreys_power,
                       method) - beta_vector_new_inner;

      beta_vector_new = beta_vector_new_inner;
      beta_vector_inner = beta_vector_new_inner;
      const double criterion_candidate =
        arma::norm(stepsize_vector_inner, "inf");
      if (criterion_inner > criterion_candidate) {
        break;
      }
    }
    criterion_vector(i - 1) = arma::norm(stepsize_vector_inner, "inf");
    beta_vector = beta_vector_new;
    eta_vector = model_matrix * beta_vector + offset;
    if (!valideta(link, arma2vec(eta_vector))) {
      Rcpp::stop(
        "invalid linear predictor - please another set of initial values for beta!!"
      );
    }
    mu_vector = linkinv(lc, eta_vector);
    if (!validmu(family, arma2vec(mu_vector))) {
      Rcpp::stop(
        "invalid fitted values - please another set of initial values for beta!!"
      );
    }
    beta_hat_matrix = arma::join_rows(beta_hat_matrix, beta_vector);
    pearson_residuals_vector =
      get_pearson_residuals(fc,
                            y_vector,
                            mu_vector,
                            weights_vector);
    if (phi_fixed == 0) {
      phi = get_phi_hat(pearson_residuals_vector, use_params);
    }
    if (alpha_fixed == 0) {
      alpha_vector = get_alpha_hat(correlation_structure,
                                   pearson_residuals_vector,
                                   id_vector,
                                   repeated_vector,
                                   phi,
                                   use_params,
                                   mdependence);
    }
    if (criterion_vector(i - 1) <= tolerance) {
      break;
    }
  }
  Rcpp::List cov_matrices =
    get_covariance_matrices_cc(y_vector,
                               model_matrix,
                               id_vector,
                               repeated_vector,
                               weights_vector,
                               link,
                               family,
                               mu_vector,
                               eta_vector,
                               correlation_structure,
                               alpha_vector,
                               phi);
  Rcpp::List ans;
  ans["beta_hat"] = beta_vector;
  ans["beta_mat"] = beta_hat_matrix;
  ans["alpha"] = alpha_vector;
  ans["phi"] = phi;
  ans["naive_covariance"] = cov_matrices[0];
  ans["robust_covariance"] = cov_matrices[1];
  ans["bc_covariance"] = cov_matrices[2];
  ans["criterion"] = criterion_vector;
  ans["eta"] = eta_vector;
  ans["residuals"] = y_vector - mu_vector;
  ans["fitted"] = mu_vector;
  ans["offset"] = offset;
  return ans;
}
//==============================================================================
