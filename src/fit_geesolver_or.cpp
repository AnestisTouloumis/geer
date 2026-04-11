#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link_functions.h"
#include "link_utils.h"
#include "utils.h"
#include "nuisance_quantities_or.h"
#include "covariance_matrices.h"
#include "clusterutils.h"
#include "method_dispatch.h"
#include <cmath>


//============================ update beta - gee OR ============================
arma::vec update_beta_gee_or(const arma::vec& y_vector,
                             const arma::mat& model_matrix,
                             const arma::vec& id_vector,
                             const arma::vec& repeated_vector,
                             const arma::vec& weights_vector,
                             const char* link,
                             const arma::vec& beta_vector,
                             const arma::vec& mu_vector,
                             const arma::vec& eta_vector,
                             const arma::vec& alpha_vector) {
  const LinkCode lc = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max = static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat naive_matrix_inverse(params_no, params_no, arma::fill::zeros);
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
    if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
      d_matrix_i.set_size(m, params_no);
    }
    d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector.subvec(a, b);
    const arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector.subvec(a, b),
                                       repeated_max,
                                       alpha_vector);
    v_matrix_i = get_v_matrix_or(mu_vector.subvec(a, b),
                                 odds_ratios_vector_i,
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
      Rcpp::stop("update_beta_gee_or: failed to solve V_i^{-1} D_i.");
    }
    if (d_matrix_trans_v_matrix_inverse_i.n_rows != params_no ||
        d_matrix_trans_v_matrix_inverse_i.n_cols != m) {
      d_matrix_trans_v_matrix_inverse_i.set_size(params_no, m);
    }
    d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector.subvec(a, b);
  }
  symmetrize_if_close(naive_matrix_inverse, 1e-10);
  arma::vec update_step(params_no, arma::fill::zeros);
  const bool ok_step = arma::solve(update_step,
                                   naive_matrix_inverse,
                                   u_vector,
                                   arma::solve_opts::likely_sympd);
  if (!ok_step) {
    Rcpp::stop("update_beta_gee_or: failed to solve for the beta update.");
  }
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - naive OR ==========================
arma::vec update_beta_naive_or(const arma::vec& y_vector,
                               const arma::mat& model_matrix,
                               const arma::vec& id_vector,
                               const arma::vec& repeated_vector,
                               const arma::vec& weights_vector,
                               const char* link,
                               const arma::vec& beta_vector,
                               const arma::vec& mu_vector,
                               const arma::vec& eta_vector,
                               const arma::vec& alpha_vector) {
  const LinkCode lc = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max = static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat lambda_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat naive_matrix_inverse(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = mueta(lc, eta_vector);
  const arma::vec delta_star_vector =
    mueta2(lc, eta_vector) / arma::square(delta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::mat d_matrix_i;
  arma::mat d_matrix_trans_i;
  arma::mat v_matrix_inverse_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::mat weights_matrix_sq_inverse_i;
  arma::mat v_matrix_tilde_inverse_i;
  arma::mat g_matrix_i;
  arma::mat identity_matrix_i;
  arma::mat kappa_matrix_delta_star_matrix_i;
  arma::vec v_matrix_tilde_inverse_mu_vector_i;
  arma::mat h_epsilon_trans_i;
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::uword m = cl.end - cl.start;
    const arma::vec mu_vector_i = mu_vector.subvec(a, b);
    const arma::vec weights_vector_i = weights_vector.subvec(a, b);
    const arma::vec s_vector_i = s_vector.subvec(a, b);
    if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
      d_matrix_i.set_size(m, params_no);
    }
    d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector.subvec(a, b);
    if (d_matrix_trans_i.n_rows != params_no || d_matrix_trans_i.n_cols != m) {
      d_matrix_trans_i.set_size(params_no, m);
    }
    d_matrix_trans_i = d_matrix_i.t();
    const arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector.subvec(a, b),
                                       repeated_max,
                                       alpha_vector);
    v_matrix_inverse_i =
      get_v_matrix_inverse_or(mu_vector_i,
                              odds_ratios_vector_i,
                              weights_vector_i);
    if (d_matrix_trans_v_matrix_inverse_i.n_rows != params_no ||
        d_matrix_trans_v_matrix_inverse_i.n_cols != m) {
      d_matrix_trans_v_matrix_inverse_i.set_size(params_no, m);
    }
    d_matrix_trans_v_matrix_inverse_i = d_matrix_trans_i * v_matrix_inverse_i;
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    const arma::vec u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector_i;
    u_vector += u_vector_i;
    weights_matrix_sq_inverse_i = arma::diagmat(1.0 / arma::sqrt(weights_vector_i));
    v_matrix_tilde_inverse_i = v_matrix_inverse_i * weights_matrix_sq_inverse_i;
    g_matrix_i = get_g_matrix(mu_vector_i, odds_ratios_vector_i);
    identity_matrix_i = arma::eye(m, m);
    kappa_matrix_delta_star_matrix_i =
      kappa_right(arma::diagmat(delta_star_vector.subvec(a, b)));
    v_matrix_tilde_inverse_mu_vector_i =
      v_matrix_tilde_inverse_i * mu_vector_i;
    h_epsilon_trans_i =
      kappa_matrix_delta_star_matrix_i +
      (arma::vectorise(v_matrix_tilde_inverse_i.t()) * mu_vector_i.t() -
      kronecker_left_identity_kappa(v_matrix_tilde_inverse_i) * g_matrix_i -
      kronecker_left_identity_kappa(
        v_matrix_tilde_inverse_i * (g_matrix_i.t() + identity_matrix_i)
      ) +
        arma::kron(v_matrix_tilde_inverse_mu_vector_i, identity_matrix_i)) *
        weights_matrix_sq_inverse_i;
    lambda_matrix -=
      arma::kron(d_matrix_trans_i, d_matrix_trans_i) *
      (arma::kron(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i.t()) -
      kronecker_sum_same(v_matrix_inverse_i) * kappa_matrix_delta_star_matrix_i -
      arma::kron(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i) +
      h_epsilon_trans_i * v_matrix_inverse_i +
      arma::kron(v_matrix_inverse_i, v_matrix_inverse_i) *
      get_v_matrix_mu_or(mu_vector_i,
                         odds_ratios_vector_i,
                         weights_vector_i)) *
                           d_matrix_i;
  }
  symmetrize_if_close(naive_matrix_inverse, 1e-10);
  arma::vec lambda_vector(params_no, arma::fill::zeros);
  for (arma::uword r = 0; r < params_no; ++r) {
    const arma::mat block =
      lambda_matrix.rows(r * params_no, (r + 1) * params_no - 1);
    const arma::mat solved_block =
      solve_chol_or_lu_mat(naive_matrix_inverse, block);
    lambda_vector[r] = 0.5 * arma::trace(solved_block);
  }
  const arma::vec update_step =
    solve_chol_or_lu_vec(naive_matrix_inverse, u_vector + lambda_vector);
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - robust OR =========================
arma::vec update_beta_robust_or(const arma::vec& y_vector,
                                const arma::mat& model_matrix,
                                const arma::vec& id_vector,
                                const arma::vec& repeated_vector,
                                const arma::vec& weights_vector,
                                const char* link,
                                const arma::vec& beta_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& eta_vector,
                                const arma::vec& alpha_vector) {
  const LinkCode lc = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max = static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat partial_derivatives_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat second_derivatives_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat naive_matrix_inverse(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = mueta(lc, eta_vector);
  const arma::vec delta_star_vector =
    mueta2(lc, eta_vector) / arma::square(delta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::mat d_matrix_i;
  arma::mat d_matrix_trans_i;
  arma::mat v_matrix_inverse_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::mat weights_matrix_sq_inverse_i;
  arma::mat v_matrix_tilde_inverse_i;
  arma::mat g_matrix_i;
  arma::mat identity_matrix_i;
  arma::mat kappa_matrix_delta_star_matrix_i;
  arma::vec v_matrix_tilde_inverse_mu_vector_i;
  arma::mat kron_d_matrix_trans_d_matrix_trans_i;
  arma::mat h_epsilon_trans_i;
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::uword m = cl.end - cl.start;
    const arma::vec mu_vector_i = mu_vector.subvec(a, b);
    const arma::vec weights_vector_i = weights_vector.subvec(a, b);
    const arma::vec s_vector_i = s_vector.subvec(a, b);
    if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
      d_matrix_i.set_size(m, params_no);
    }
    d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector.subvec(a, b);
    if (d_matrix_trans_i.n_rows != params_no || d_matrix_trans_i.n_cols != m) {
      d_matrix_trans_i.set_size(params_no, m);
    }
    d_matrix_trans_i = d_matrix_i.t();
    const arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector.subvec(a, b),
                                       repeated_max,
                                       alpha_vector);
    v_matrix_inverse_i =
      get_v_matrix_inverse_or(mu_vector_i,
                              odds_ratios_vector_i,
                              weights_vector_i);
    if (d_matrix_trans_v_matrix_inverse_i.n_rows != params_no ||
        d_matrix_trans_v_matrix_inverse_i.n_cols != m) {
      d_matrix_trans_v_matrix_inverse_i.set_size(params_no, m);
    }
    d_matrix_trans_v_matrix_inverse_i = d_matrix_trans_i * v_matrix_inverse_i;
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    const arma::vec u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector_i;
    u_vector += u_vector_i;
    meat_matrix += u_vector_i * u_vector_i.t();
    weights_matrix_sq_inverse_i = arma::diagmat(1.0 / arma::sqrt(weights_vector_i));
    v_matrix_tilde_inverse_i = v_matrix_inverse_i * weights_matrix_sq_inverse_i;
    g_matrix_i = get_g_matrix(mu_vector_i, odds_ratios_vector_i);
    identity_matrix_i = arma::eye(m, m);
    kappa_matrix_delta_star_matrix_i =
      kappa_right(arma::diagmat(delta_star_vector.subvec(a, b)));
    v_matrix_tilde_inverse_mu_vector_i =
      v_matrix_tilde_inverse_i * mu_vector_i;
    kron_d_matrix_trans_d_matrix_trans_i =
      arma::kron(d_matrix_trans_i, d_matrix_trans_i);
    h_epsilon_trans_i =
      kappa_matrix_delta_star_matrix_i +
      (arma::vectorise(v_matrix_tilde_inverse_i.t()) * mu_vector_i.t() -
      kronecker_left_identity_kappa(v_matrix_tilde_inverse_i) * g_matrix_i -
      kronecker_left_identity_kappa(
        v_matrix_tilde_inverse_i * (g_matrix_i.t() + identity_matrix_i)
      ) +
        arma::kron(v_matrix_tilde_inverse_mu_vector_i, identity_matrix_i)) *
        weights_matrix_sq_inverse_i;
    partial_derivatives_matrix +=
      kron_d_matrix_trans_d_matrix_trans_i *
      h_epsilon_trans_i *
      v_matrix_inverse_i *
      s_vector_i *
      u_vector_i.t();
    second_derivatives_matrix +=
      kron_d_matrix_trans_d_matrix_trans_i *
      (arma::kron(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i.t()) -
      kronecker_sum_same(v_matrix_inverse_i) * kappa_matrix_delta_star_matrix_i -
      arma::kron(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i) -
      h_epsilon_trans_i * v_matrix_inverse_i +
      arma::kron(v_matrix_inverse_i, v_matrix_inverse_i) *
      get_v_matrix_mu_or(mu_vector_i,
                         odds_ratios_vector_i,
                         weights_vector_i)) *
                           d_matrix_i;
  }
  symmetrize_if_close(naive_matrix_inverse, 1e-10);
  const arma::mat naive_matrix_meat_matrix =
    solve_chol_or_lu_mat(naive_matrix_inverse, meat_matrix);
  const arma::mat robust_matrix =
    solve_chol_or_lu_mat(naive_matrix_inverse, naive_matrix_meat_matrix.t());
  arma::vec lambda_vector(params_no, arma::fill::zeros);
  for (arma::uword r = 0; r < params_no; ++r) {
    const arma::mat first_block =
      partial_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
    const arma::mat second_block =
      second_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
    const arma::mat solved_first_block =
      solve_chol_or_lu_mat(naive_matrix_inverse, first_block);
    lambda_vector[r] =
      -(arma::trace(solved_first_block) +
      0.5 * arma::trace(robust_matrix * second_block));
  }
  const arma::vec update_step =
    solve_chol_or_lu_vec(naive_matrix_inverse, u_vector + lambda_vector);
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - empirical OR ======================
arma::vec update_beta_empirical_or(const arma::vec& y_vector,
                                   const arma::mat& model_matrix,
                                   const arma::vec& id_vector,
                                   const arma::vec& repeated_vector,
                                   const arma::vec& weights_vector,
                                   const char* link,
                                   const arma::vec& beta_vector,
                                   const arma::vec& mu_vector,
                                   const arma::vec& eta_vector,
                                   const arma::vec& alpha_vector) {
  const LinkCode lc = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max = static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat partial_derivatives_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat second_derivatives_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat observed_fisher_info_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = mueta(lc, eta_vector);
  const arma::vec mueta2_vector = mueta2(lc, eta_vector);
  const arma::vec delta_star_vector =
    mueta2_vector / arma::square(delta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const arma::vec delta_tilde_star_vector =
    (delta_vector % mueta3(lc, eta_vector) - 2.0 * arma::square(mueta2_vector)) /
      arma::pow(delta_vector, 4.0);
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::mat d_matrix_i;
  arma::mat d_matrix_trans_i;
  arma::mat v_matrix_inverse_i;
  arma::mat delta_star_matrix_i;
  arma::mat weights_matrix_sq_inverse_i;
  arma::mat w_matrix_i;
  arma::mat identity_matrix_i;
  arma::mat v_matrix_inverse_derivative_i;
  arma::mat epsilon_matrix_i;
  arma::mat observed_fisher_info_matrix_i;
  arma::mat v_matrix_tilde_inverse_i;
  arma::mat g_matrix_i;
  arma::mat h_epsilon_matrix_trans_i;
  arma::mat epsilon_matrix_transpose_derivative_term1_i;
  arma::mat epsilon_matrix_transpose_derivative_term2_i;
  arma::mat w_tilde_matrix_i;
  arma::mat epsilon_matrix_transpose_derivative_term3_i;
  arma::mat epsilon_matrix_transpose_derivative_term4_i;
  arma::mat epsilon_matrix_transpose_derivative_i;
  arma::mat second_derivatives_matrix_terms12_i;
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::uword m = cl.end - cl.start;
    const arma::vec mu_vector_i = mu_vector.subvec(a, b);
    const arma::vec s_vector_i = s_vector.subvec(a, b);
    const arma::vec weights_vector_i = weights_vector.subvec(a, b);
    if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
      d_matrix_i.set_size(m, params_no);
    }
    d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector.subvec(a, b);
    if (d_matrix_trans_i.n_rows != params_no || d_matrix_trans_i.n_cols != m) {
      d_matrix_trans_i.set_size(params_no, m);
    }
    d_matrix_trans_i = d_matrix_i.t();
    const arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector.subvec(a, b),
                                       repeated_max,
                                       alpha_vector);
    v_matrix_inverse_i =
      get_v_matrix_inverse_or(mu_vector_i,
                              odds_ratios_vector_i,
                              weights_vector_i);
    const arma::vec u_vector_i = d_matrix_trans_i * v_matrix_inverse_i * s_vector_i;
    u_vector += u_vector_i;
    meat_matrix += u_vector_i * u_vector_i.t();
    delta_star_matrix_i = arma::diagmat(delta_star_vector.subvec(a, b));
    weights_matrix_sq_inverse_i = arma::diagmat(1.0 / arma::sqrt(weights_vector_i));
    w_matrix_i = arma::diagmat(v_matrix_inverse_i * s_vector_i);
    identity_matrix_i = arma::eye(m, m);
    v_matrix_inverse_derivative_i =
      -arma::kron(v_matrix_inverse_i, v_matrix_inverse_i) *
      get_v_matrix_mu_or(mu_vector_i,
                         odds_ratios_vector_i,
                         weights_vector_i);
    epsilon_matrix_i =
      delta_star_matrix_i * w_matrix_i -
      v_matrix_inverse_i +
      arma::kron(s_vector_i.t(), identity_matrix_i) * v_matrix_inverse_derivative_i;
    observed_fisher_info_matrix_i =
      d_matrix_trans_i * epsilon_matrix_i * d_matrix_i;
    observed_fisher_info_matrix -= observed_fisher_info_matrix_i;
    partial_derivatives_matrix +=
      arma::vectorise(observed_fisher_info_matrix_i.t()) * u_vector_i.t();
    v_matrix_tilde_inverse_i =
      v_matrix_inverse_i * weights_matrix_sq_inverse_i;
    const arma::vec v_matrix_tilde_inverse_trans_s_vector_i =
      v_matrix_tilde_inverse_i.t() * s_vector_i;
    const arma::vec v_matrix_tilde_inverse_mu_vector_i =
      v_matrix_tilde_inverse_i * mu_vector_i;
    epsilon_matrix_transpose_derivative_term1_i =
      kappa_right(w_matrix_i * arma::diagmat(delta_tilde_star_vector.subvec(a, b))) +
      arma::vectorise(v_matrix_tilde_inverse_i.t()) * v_matrix_tilde_inverse_trans_s_vector_i.t() -
      kronecker_vector_matrix(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i) +
      arma::kron(v_matrix_tilde_inverse_i, v_matrix_tilde_inverse_trans_s_vector_i);
    g_matrix_i = get_g_matrix(mu_vector_i, odds_ratios_vector_i);
    epsilon_matrix_transpose_derivative_term2_i =
      kappa_right(delta_star_matrix_i) +
      (arma::vectorise(v_matrix_tilde_inverse_i.t()) * mu_vector_i.t() -
      kronecker_left_identity_kappa(v_matrix_tilde_inverse_i) * g_matrix_i -
      kronecker_left_identity_kappa(
        v_matrix_tilde_inverse_i * (g_matrix_i.t() + identity_matrix_i)
      )) *
        weights_matrix_sq_inverse_i;
    epsilon_matrix_transpose_derivative_term2_i =
      epsilon_matrix_transpose_derivative_term2_i *
      (epsilon_matrix_i - w_matrix_i * delta_star_matrix_i);
    w_tilde_matrix_i = weights_matrix_sq_inverse_i * w_matrix_i;
    h_epsilon_matrix_trans_i =
      v_matrix_tilde_inverse_trans_s_vector_i * mu_vector_i.t() -
      w_tilde_matrix_i * (identity_matrix_i + g_matrix_i) +
      arma::diagmat(
        -g_matrix_i * v_matrix_tilde_inverse_trans_s_vector_i +
          arma::as_scalar(v_matrix_tilde_inverse_trans_s_vector_i.t() * mu_vector_i)
      );
    epsilon_matrix_transpose_derivative_term3_i =
      arma::kron(v_matrix_tilde_inverse_mu_vector_i * s_vector_i.t(),
                 weights_matrix_sq_inverse_i) +
                   arma::kron(identity_matrix_i,
                              weights_matrix_sq_inverse_i * h_epsilon_matrix_trans_i - identity_matrix_i);
    epsilon_matrix_transpose_derivative_term3_i =
      epsilon_matrix_transpose_derivative_term3_i * v_matrix_inverse_derivative_i;
    epsilon_matrix_transpose_derivative_term4_i =
      (arma::kron(v_matrix_tilde_inverse_i, w_tilde_matrix_i) +
      kronecker_left_identity_kappa(v_matrix_tilde_inverse_i) *
      kronecker_vector_identity(v_matrix_tilde_inverse_trans_s_vector_i).t()) *
      get_g_matrix_mu(mu_vector_i, odds_ratios_vector_i);
    epsilon_matrix_transpose_derivative_i =
      epsilon_matrix_transpose_derivative_term1_i +
      epsilon_matrix_transpose_derivative_term2_i +
      epsilon_matrix_transpose_derivative_term3_i -
      epsilon_matrix_transpose_derivative_term4_i;
    second_derivatives_matrix_terms12_i =
      kronecker_left_identity_kappa(epsilon_matrix_i * delta_star_matrix_i) +
      kronecker_identity_right_kappa(epsilon_matrix_i.t() * delta_star_matrix_i);
    second_derivatives_matrix +=
      arma::kron(d_matrix_trans_i, d_matrix_trans_i) *
      (second_derivatives_matrix_terms12_i + epsilon_matrix_transpose_derivative_i) *
      d_matrix_i;
  }
  arma::mat robust_matrix(params_no, params_no, arma::fill::zeros);
  arma::vec lambda_vector(params_no, arma::fill::zeros);
  arma::mat lu_lower, lu_upper, permutation_matrix;
  const bool lu_success =
    arma::lu(lu_lower, lu_upper, permutation_matrix, observed_fisher_info_matrix);
  if (!lu_success) {
    Rcpp::stop("update_beta_empirical_or: LU factorization failed.");
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
      Rcpp::stop("update_beta_empirical_or: forward LU solve failed.");
    }
    arma::mat ans;
    const bool ok_back =
      arma::solve(ans,
                  arma::trimatu(lu_upper),
                  lu_forward,
                  arma::solve_opts::allow_ugly);
    if (!ok_back) {
      Rcpp::stop("update_beta_empirical_or: backward LU solve failed.");
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
      Rcpp::stop("update_beta_empirical_or: forward LU solve failed.");
    }
    arma::vec ans;
    const bool ok_back =
      arma::solve(ans,
                  arma::trimatu(lu_upper),
                  lu_forward,
                  arma::solve_opts::allow_ugly);
    if (!ok_back) {
      Rcpp::stop("update_beta_empirical_or: backward LU solve failed.");
    }
    return ans;
  };
  const arma::mat fisher_inv_meat = solve_observed_fisher_matrix(meat_matrix);
  robust_matrix = solve_observed_fisher_matrix(fisher_inv_meat.t());
  for (arma::uword r = 0; r < params_no; ++r) {
    const arma::mat first_block =
      partial_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
    const arma::mat second_block =
      second_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
    const arma::mat solved_first_block =
      solve_observed_fisher_matrix(first_block);
    lambda_vector[r] =
      -(arma::trace(solved_first_block) +
      0.5 * arma::trace(robust_matrix * second_block));
  }
  const arma::vec update_step =
    solve_observed_fisher_vector(u_vector + lambda_vector);
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - jeffreys OR =======================
arma::vec update_beta_jeffreys_or(const arma::vec& y_vector,
                                  const arma::mat& model_matrix,
                                  const arma::vec& id_vector,
                                  const arma::vec& repeated_vector,
                                  const arma::vec& weights_vector,
                                  const char* link,
                                  const arma::vec& beta_vector,
                                  const arma::vec& mu_vector,
                                  const arma::vec& eta_vector,
                                  const arma::vec& alpha_vector,
                                  const double& jeffreys_power) {
  const LinkCode lc = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max = static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat naive_matrix_inverse(params_no, params_no, arma::fill::zeros);
  arma::mat naive_matrix_inverse_derivative(params_no * params_no,
                                            params_no,
                                            arma::fill::zeros);
  const arma::vec delta_vector = mueta(lc, eta_vector);
  const arma::vec delta_star_vector =
    mueta2(lc, eta_vector) / arma::square(delta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const auto clusters = clusters_from_sorted_id(id_vector);
  arma::mat d_matrix_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::mat v_matrix_inverse_i;
  for (const auto& cl : clusters) {
    const arma::uword a = cl.start;
    const arma::uword b = cl.end - 1;
    const arma::uword m = cl.end - cl.start;
    const arma::vec mu_vector_i = mu_vector.subvec(a, b);
    const arma::vec weights_vector_i = weights_vector.subvec(a, b);
    if (d_matrix_i.n_rows != m || d_matrix_i.n_cols != params_no) {
      d_matrix_i.set_size(m, params_no);
    }
    d_matrix_i = model_matrix.rows(a, b);
    d_matrix_i.each_col() %= delta_vector.subvec(a, b);
    const arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector.subvec(a, b),
                                       repeated_max,
                                       alpha_vector);
    v_matrix_inverse_i =
      get_v_matrix_inverse_or(mu_vector_i,
                              odds_ratios_vector_i,
                              weights_vector_i);
    if (d_matrix_trans_v_matrix_inverse_i.n_rows != params_no ||
        d_matrix_trans_v_matrix_inverse_i.n_cols != m) {
      d_matrix_trans_v_matrix_inverse_i.set_size(params_no, m);
    }
    d_matrix_trans_v_matrix_inverse_i = d_matrix_i.t() * v_matrix_inverse_i;
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector.subvec(a, b);
    naive_matrix_inverse_derivative +=
      arma::kron(d_matrix_i, d_matrix_i).t() *
      (kronecker_sum_same(
          v_matrix_inverse_i * arma::diagmat(delta_star_vector.subvec(a, b))
      ) *
        kappa_matrix(m) -
        arma::kron(v_matrix_inverse_i, v_matrix_inverse_i) *
        get_v_matrix_mu_or(mu_vector_i,
                           odds_ratios_vector_i,
                           weights_vector_i)) *
                             d_matrix_i;
  }
  symmetrize_if_close(naive_matrix_inverse, 1e-10);
  const arma::mat naive_matrix =
    solve_chol_or_lu_mat(naive_matrix_inverse, arma::eye(params_no, params_no));
  const arma::vec lambda_vector =
    jeffreys_power *
    naive_matrix_inverse_derivative.t() *
    arma::vectorise(naive_matrix);
  const arma::vec update_step =
    solve_chol_or_lu_vec(naive_matrix_inverse, u_vector + lambda_vector);
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - or ================================
arma::vec update_beta_or(const arma::vec& y_vector,
                         const arma::mat& model_matrix,
                         const arma::vec& id_vector,
                         const arma::vec& repeated_vector,
                         const arma::vec& weights_vector,
                         const char* link,
                         const arma::vec& beta_vector,
                         const arma::vec& mu_vector,
                         const arma::vec& eta_vector,
                         const arma::vec& alpha_vector,
                         const double& jeffreys_power,
                         const char* method) {
  switch (method_code(method)) {
  case M_GEE:
    return update_beta_gee_or(y_vector,
                              model_matrix,
                              id_vector,
                              repeated_vector,
                              weights_vector,
                              link,
                              beta_vector,
                              mu_vector,
                              eta_vector,
                              alpha_vector);
  case M_BR_NAIVE:
    return update_beta_naive_or(y_vector,
                                model_matrix,
                                id_vector,
                                repeated_vector,
                                weights_vector,
                                link,
                                beta_vector,
                                mu_vector,
                                eta_vector,
                                alpha_vector);
  case M_BR_ROBUST:
    return update_beta_robust_or(y_vector,
                                 model_matrix,
                                 id_vector,
                                 repeated_vector,
                                 weights_vector,
                                 link,
                                 beta_vector,
                                 mu_vector,
                                 eta_vector,
                                 alpha_vector);
  case M_BR_EMPIRICAL:
    return update_beta_empirical_or(y_vector,
                                    model_matrix,
                                    id_vector,
                                    repeated_vector,
                                    weights_vector,
                                    link,
                                    beta_vector,
                                    mu_vector,
                                    eta_vector,
                                    alpha_vector);
  case M_PGEE_JEFFREYS:
    return update_beta_jeffreys_or(y_vector,
                                   model_matrix,
                                   id_vector,
                                   repeated_vector,
                                   weights_vector,
                                   link,
                                   beta_vector,
                                   mu_vector,
                                   eta_vector,
                                   alpha_vector,
                                   jeffreys_power);
  }
  Rcpp::stop("update_beta_or: unsupported method.");
}
//==============================================================================


//=========================== fitting function =================================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List fit_bingee_or(const arma::vec& y_vector,
                         const arma::mat& model_matrix,
                         const arma::vec& id_vector,
                         const arma::vec& repeated_vector,
                         const arma::vec& weights_vector,
                         const char* link,
                         arma::vec beta_vector,
                         const arma::vec& offset,
                         const int& maxiter,
                         const double& tolerance,
                         const int& step_maxiter,
                         const int& step_multiplier,
                         const double& jeffreys_power,
                         const char* method,
                         const arma::vec& alpha_vector) {
  const LinkCode lc = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
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
      "invalid linear predictor - please another set of initial values for beta!!"
    );
  }
  arma::vec mu_vector = linkinv(lc, eta_vector);
  if (!validmu("binomial", arma2vec(mu_vector))) {
    Rcpp::stop(
      "invalid fitted values - please another set of initial values for beta!!"
    );
  }
  for (int i = 1; i < maxiter + 1; ++i) {
    stepsize_vector =
      update_beta_or(y_vector,
                     model_matrix,
                     id_vector,
                     repeated_vector,
                     weights_vector,
                     link,
                     beta_vector,
                     mu_vector,
                     eta_vector,
                     alpha_vector,
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
          "invalid initial linear predictor: please another set of initial values for beta!!"
        );
      }
      mu_vector = linkinv(lc, eta_vector);
      if (!validmu("binomial", arma2vec(mu_vector))) {
        Rcpp::stop(
          "invalid fitted values - please another set of initial values for beta!!"
        );
      }
      stepsize_vector_inner =
        update_beta_or(y_vector,
                       model_matrix,
                       id_vector,
                       repeated_vector,
                       weights_vector,
                       link,
                       beta_vector_new_inner,
                       mu_vector,
                       eta_vector,
                       alpha_vector,
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
    if (!validmu("binomial", arma2vec(mu_vector))) {
      Rcpp::stop(
        "invalid fitted values - please another set of initial values for beta!!"
      );
    }
    beta_hat_matrix = arma::join_rows(beta_hat_matrix, beta_vector);
    if (criterion_vector(i - 1) <= tolerance) {
      break;
    }
  }
  Rcpp::List cov_matrices =
    get_covariance_matrices_or(y_vector,
                               model_matrix,
                               id_vector,
                               repeated_vector,
                               weights_vector,
                               link,
                               mu_vector,
                               eta_vector,
                               alpha_vector);
  Rcpp::List ans;
  ans["beta_hat"] = beta_vector;
  ans["beta_mat"] = beta_hat_matrix;
  ans["alpha"] = alpha_vector;
  ans["phi"] = 1.0;
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
