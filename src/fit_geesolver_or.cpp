#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link_functions.h"
#include "utils.h"
#include "nuisance_quantities_or.h"
#include "covariance_matrices.h"

//============================ update beta - gee OR ============================
// [[Rcpp::export]]
arma::vec update_beta_gee_or(const arma::vec & y_vector,
                             const arma::mat & model_matrix,
                             const arma::vec & id_vector,
                             const arma::vec & repeated_vector,
                             const arma::vec & weights_vector,
                             const char * link,
                             const arma::vec & beta_vector,
                             const arma::vec & mu_vector,
                             const arma::vec & eta_vector,
                             const arma::vec & alpha_vector) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  int repeated_max = max(repeated_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec s_vector = y_vector - mu_vector;
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::mat d_matrix_i =
      arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max,
                                       alpha_vector);
    arma::mat d_matrix_trans_v_matrix_inverse_i =
      trans(
        solve(get_v_matrix_or(mu_vector(id_vector_i),
                              odds_ratios_vector_i,
                              weights_vector(id_vector_i)),
                              d_matrix_i)
      );
    naive_matrix_inverse +=
      d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    u_vector +=
      d_matrix_trans_v_matrix_inverse_i * s_vector(id_vector_i);
  }
  arma::vec ans = beta_vector + solve(naive_matrix_inverse, u_vector);
  return ans;
}
//==============================================================================


//============================ update beta - naive OR ==========================
// [[Rcpp::export]]
arma::vec update_beta_naive_or(const arma::vec & y_vector,
                               const arma::mat & model_matrix,
                               const arma::vec & id_vector,
                               const arma::vec & repeated_vector,
                               const arma::vec & weights_vector,
                               const char * link,
                               const arma::vec & beta_vector,
                               const arma::vec & mu_vector,
                               const arma::vec & eta_vector,
                               const arma::vec & alpha_vector) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  int repeated_max = max(repeated_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat lambda_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec s_vector = y_vector - mu_vector;
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec mu_vector_i = mu_vector(id_vector_i);
    arma::vec weights_vector_i = weights_vector(id_vector_i);
    arma::mat d_matrix_i =
      diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::mat d_matrix_trans_i = trans(d_matrix_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max,
                                       alpha_vector);
    arma::mat v_matrix_inverse_i =
      get_v_matrix_inverse_or(mu_vector_i,
                              odds_ratios_vector_i,
                              weights_vector_i);
    arma::mat d_matrix_trans_v_matrix_inverse_i =
      d_matrix_trans_i * v_matrix_inverse_i;
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    arma::vec s_vector_i = s_vector(id_vector_i);
    arma::vec u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector_i;
    u_vector += u_vector_i;
    arma::mat weights_matrix_sq_inverse_i = arma::diagmat(1/sqrt(weights_vector_i));
    arma::mat v_matrix_tilde_inverse_i =
      v_matrix_inverse_i * weights_matrix_sq_inverse_i;
    arma::mat g_matrix_i = get_g_matrix(mu_vector_i, odds_ratios_vector_i);
    double cluster_size_i = id_vector_i.n_elem;
    arma::mat identity_matrix_i = arma::eye(cluster_size_i, cluster_size_i);
    arma::mat kappa_matrix_delta_star_matrix_i =
      kappa_right(arma::diagmat(delta_star_vector(id_vector_i)));
    arma::mat v_matrix_tilde_inverse_mu_vector_i =
      v_matrix_tilde_inverse_i * mu_vector_i;
    arma::mat h_epsilon_trans_i =
      kappa_matrix_delta_star_matrix_i +
      (
          arma::vectorise(trans(v_matrix_tilde_inverse_i)) * trans(mu_vector_i) -
            kronecker_left_identity_kappa(v_matrix_tilde_inverse_i) * g_matrix_i -
            kronecker_left_identity_kappa(v_matrix_tilde_inverse_i * (trans(g_matrix_i) + identity_matrix_i)) +
            kron(v_matrix_tilde_inverse_mu_vector_i, identity_matrix_i)
      ) * weights_matrix_sq_inverse_i;
    lambda_matrix -=
      kron(d_matrix_trans_i, d_matrix_trans_i) *
      (
          kron(v_matrix_tilde_inverse_mu_vector_i, trans(v_matrix_tilde_inverse_i)) -
            kronecker_sum_same(v_matrix_inverse_i) * kappa_matrix_delta_star_matrix_i -
            kron(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i) +
            h_epsilon_trans_i * v_matrix_inverse_i +
            kron(v_matrix_inverse_i, v_matrix_inverse_i) *
            get_v_matrix_mu_or(mu_vector_i,
                               odds_ratios_vector_i,
                               weights_vector_i)
      ) *
        d_matrix_i;
  }
  arma::vec lambda_vector = arma::zeros(params_no);
  for(int r = 1; r < params_no + 1; r++) {
    lambda_vector(r - 1) =
      0.5 * trace(solve(naive_matrix_inverse,
                        lambda_matrix.rows((r - 1) * params_no, r * params_no - 1)));
  }
  arma::vec ans = beta_vector +
    solve(naive_matrix_inverse, u_vector + lambda_vector);
  return ans;
}
//==============================================================================


//============================ update beta - robust OR =========================
// [[Rcpp::export]]
arma::vec update_beta_robust_or(const arma::vec & y_vector,
                                const arma::mat & model_matrix,
                                const arma::vec & id_vector,
                                const arma::vec & repeated_vector,
                                const arma::vec & weights_vector,
                                const char * link,
                                const arma::vec & beta_vector,
                                const arma::vec & mu_vector,
                                const arma::vec & eta_vector,
                                const arma::vec & alpha_vector) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  int repeated_max = max(repeated_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat partial_derivatives_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat second_derivatives_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::mat meat_matrix = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec s_vector = y_vector - mu_vector;
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec mu_vector_i = mu_vector(id_vector_i);
    arma::vec weights_vector_i = weights_vector(id_vector_i);
    arma::mat d_matrix_i =
      diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::mat d_matrix_trans_i = trans(d_matrix_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max,
                                       alpha_vector);
    arma::mat v_matrix_inverse_i =
      get_v_matrix_inverse_or(mu_vector_i,
                              odds_ratios_vector_i,
                              weights_vector_i);
    arma::mat d_matrix_trans_v_matrix_inverse_i =
      d_matrix_trans_i * v_matrix_inverse_i;
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    arma::vec s_vector_i = s_vector(id_vector_i);
    arma::vec u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector_i;
    u_vector += u_vector_i;
    meat_matrix += u_vector_i * trans(u_vector_i);
    arma::mat weights_matrix_sq_inverse_i = arma::diagmat(1/sqrt(weights_vector_i));
    arma::mat v_matrix_tilde_inverse_i =
      v_matrix_inverse_i * weights_matrix_sq_inverse_i;
    arma::mat g_matrix_i = get_g_matrix(mu_vector_i, odds_ratios_vector_i);
    double cluster_size_i = id_vector_i.n_elem;
    arma::mat identity_matrix_i = arma::eye(cluster_size_i, cluster_size_i);
    arma::mat kappa_matrix_delta_star_matrix_i =
      kappa_right(arma::diagmat(delta_star_vector(id_vector_i)));
    arma::mat v_matrix_tilde_inverse_mu_vector_i =
      v_matrix_tilde_inverse_i * mu_vector_i;
    arma::mat kron_d_matrix_trans_d_matrix_trans_i =
      kron(d_matrix_trans_i, d_matrix_trans_i);
    arma::mat h_epsilon_trans_i =
      kappa_matrix_delta_star_matrix_i +
      (
          arma::vectorise(trans(v_matrix_tilde_inverse_i)) * trans(mu_vector_i) -
            kronecker_left_identity_kappa(v_matrix_tilde_inverse_i) * g_matrix_i -
            kronecker_left_identity_kappa(v_matrix_tilde_inverse_i * (trans(g_matrix_i) + identity_matrix_i)) +
            kron(v_matrix_tilde_inverse_mu_vector_i, identity_matrix_i)
      ) * weights_matrix_sq_inverse_i;
    partial_derivatives_matrix +=
      kron_d_matrix_trans_d_matrix_trans_i *
      h_epsilon_trans_i *
      v_matrix_inverse_i *
      s_vector_i *
      trans(u_vector_i);
    second_derivatives_matrix +=
      kron_d_matrix_trans_d_matrix_trans_i *
      (
          kron(v_matrix_tilde_inverse_mu_vector_i, trans(v_matrix_tilde_inverse_i)) -
            kronecker_sum_same(v_matrix_inverse_i) * kappa_matrix_delta_star_matrix_i -
            kron(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i) -
            h_epsilon_trans_i * v_matrix_inverse_i +
            kron(v_matrix_inverse_i, v_matrix_inverse_i) *
            get_v_matrix_mu_or(mu_vector_i,
                               odds_ratios_vector_i,
                               weights_vector_i)
      ) *
        d_matrix_i;
  }
  arma::mat robust_matrix =
    solve(naive_matrix_inverse,
          trans(
            solve(naive_matrix_inverse,
                  meat_matrix)
          )
    );
  arma::vec lambda_vector = arma::zeros(params_no);
  for(int r = 1; r < params_no + 1; r++) {
    lambda_vector(r - 1) = -
      (trace(
          solve(naive_matrix_inverse,
                partial_derivatives_matrix.rows((r - 1) * params_no,
                                                r * params_no - 1)
          )
      ) +
        0.5 * trace(
            robust_matrix *
              second_derivatives_matrix.rows((r - 1) * params_no, r * params_no - 1)
        )
      );
  }
  arma::vec ans =
    beta_vector + solve(naive_matrix_inverse, u_vector + lambda_vector);
  return ans;
}
//==============================================================================


//============================ update beta - empirical OR ======================
// [[Rcpp::export]]
arma::vec update_beta_empirical_or(const arma::vec & y_vector,
                                   const arma::mat & model_matrix,
                                   const arma::vec & id_vector,
                                   const arma::vec & repeated_vector,
                                   const arma::vec & weights_vector,
                                   const char * link,
                                   const arma::vec & beta_vector,
                                   const arma::vec & mu_vector,
                                   const arma::vec & eta_vector,
                                   const arma::vec & alpha_vector) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  int repeated_max = max(repeated_vector);
  arma::vec u_vector = arma::zeros(params_no);
  arma::mat partial_derivatives_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat second_derivatives_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat observed_fisher_info_matrix = arma::zeros(params_no, params_no);
  arma::mat meat_matrix = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec s_vector = y_vector - mu_vector;
  arma::vec delta_tilde_star_vector =
    (delta_vector % mueta3(link, eta_vector) - 2 * pow(mueta2(link, eta_vector), 2))/
      pow(delta_vector, 4);
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec mu_vector_i = mu_vector(id_vector_i);
    arma::vec eta_vector_i = eta_vector(id_vector_i);
    arma::vec mueta_vector_i = mueta(link, eta_vector_i);
    arma::mat d_matrix_i =
      arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::mat d_matrix_trans_i= trans(d_matrix_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max,
                                       alpha_vector);
    arma::vec weights_vector_i = weights_vector(id_vector_i);
    arma::mat v_matrix_inverse_i =
      get_v_matrix_inverse_or(mu_vector_i,
                              odds_ratios_vector_i,
                              weights_vector_i);
    arma::vec s_vector_i = s_vector(id_vector_i);
    arma::vec u_vector_i = d_matrix_trans_i * v_matrix_inverse_i * s_vector_i;
    u_vector += u_vector_i;
    meat_matrix += u_vector_i * trans(u_vector_i);
    arma::mat delta_star_matrix_i = arma::diagmat(delta_star_vector(id_vector_i));
    arma::mat weights_matrix_sq_inverse_i = arma::diagmat(1/sqrt(weights_vector_i));
    arma::mat w_matrix_i = arma::diagmat(v_matrix_inverse_i * s_vector_i);
    double cluster_size_i = id_vector_i.n_elem;
    arma::mat identity_matrix_i = arma::eye(cluster_size_i, cluster_size_i);
    arma::mat v_matrix_inverse_derivative_i =
      - kron(v_matrix_inverse_i, v_matrix_inverse_i) *
      get_v_matrix_mu_or(mu_vector_i,
                         odds_ratios_vector_i,
                         weights_vector_i);
    arma::mat epsilon_matrix_i =
      delta_star_matrix_i * w_matrix_i -
      v_matrix_inverse_i +
      kron(trans(s_vector_i), identity_matrix_i) * v_matrix_inverse_derivative_i;
    arma::mat observed_fisher_info_matrix_i =
      d_matrix_trans_i * epsilon_matrix_i * d_matrix_i;
    observed_fisher_info_matrix -= observed_fisher_info_matrix_i;
    partial_derivatives_matrix +=
      arma::vectorise(trans(observed_fisher_info_matrix_i)) * trans(u_vector_i);
    // first row
    arma::mat v_matrix_tilde_inverse_i =
      v_matrix_inverse_i * weights_matrix_sq_inverse_i;
    arma::vec v_matrix_tilde_inverse_trans_s_vector_i =
      trans(v_matrix_tilde_inverse_i) * s_vector_i;
    arma::vec v_matrix_tilde_inverse_mu_vector_i =
      v_matrix_tilde_inverse_i * mu_vector_i;
    arma::mat epsilon_matrix_transpose_derivative_term1_i =
      kappa_right(w_matrix_i * arma::diagmat(delta_tilde_star_vector(id_vector_i))) +
      arma::vectorise(trans(v_matrix_tilde_inverse_i)) * trans(v_matrix_tilde_inverse_trans_s_vector_i) -
      kronecker_vector_matrix(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i) +
      kron(v_matrix_tilde_inverse_i, v_matrix_tilde_inverse_trans_s_vector_i);
    //second row
    arma::mat g_matrix_i = get_g_matrix(mu_vector_i, odds_ratios_vector_i);
    arma::mat epsilon_matrix_transpose_derivative_term2_i =
      kappa_right(delta_star_matrix_i) +
      (
          arma::vectorise(trans(v_matrix_tilde_inverse_i)) * trans(mu_vector_i) -
            kronecker_left_identity_kappa(v_matrix_tilde_inverse_i) * g_matrix_i -
            kronecker_left_identity_kappa(v_matrix_tilde_inverse_i * (trans(g_matrix_i) + identity_matrix_i))
      ) * weights_matrix_sq_inverse_i;
    epsilon_matrix_transpose_derivative_term2_i =
      epsilon_matrix_transpose_derivative_term2_i * (epsilon_matrix_i - w_matrix_i * delta_star_matrix_i);
    // third row
    arma::mat w_tilde_matrix_i = weights_matrix_sq_inverse_i * w_matrix_i;
    arma::mat h_epsilon_matrix_trans_i =
      v_matrix_tilde_inverse_trans_s_vector_i * trans(mu_vector_i) -
      w_tilde_matrix_i * (identity_matrix_i + g_matrix_i) +
      arma::diagmat(- g_matrix_i * v_matrix_tilde_inverse_trans_s_vector_i +
      as_scalar(trans(v_matrix_tilde_inverse_trans_s_vector_i) * mu_vector_i));
    arma::mat epsilon_matrix_transpose_derivative_term3_i =
      kron(v_matrix_tilde_inverse_mu_vector_i * trans(s_vector_i),
           weights_matrix_sq_inverse_i) +
             kron(identity_matrix_i,
                  weights_matrix_sq_inverse_i * h_epsilon_matrix_trans_i - identity_matrix_i);
    epsilon_matrix_transpose_derivative_term3_i =
      epsilon_matrix_transpose_derivative_term3_i * v_matrix_inverse_derivative_i;
    // fourth row
    arma::mat epsilon_matrix_transpose_derivative_term4_i =
      (kron(v_matrix_tilde_inverse_i, w_tilde_matrix_i) +
      kronecker_left_identity_kappa(v_matrix_tilde_inverse_i) *
      trans(kronecker_vector_identity(v_matrix_tilde_inverse_trans_s_vector_i))) *
      get_g_matrix_mu(mu_vector_i, odds_ratios_vector_i);
    arma::mat epsilon_matrix_transpose_derivative_i =
      epsilon_matrix_transpose_derivative_term1_i +
      epsilon_matrix_transpose_derivative_term2_i +
      epsilon_matrix_transpose_derivative_term3_i -
      epsilon_matrix_transpose_derivative_term4_i;
    arma::mat second_derivatives_matrix_terms12_i =
      kronecker_left_identity_kappa(epsilon_matrix_i * delta_star_matrix_i) +
      kronecker_identity_right_kappa(trans(epsilon_matrix_i) * delta_star_matrix_i);
    second_derivatives_matrix +=
      kron(d_matrix_trans_i, d_matrix_trans_i) *
      (second_derivatives_matrix_terms12_i + epsilon_matrix_transpose_derivative_i) *
      d_matrix_i;
  }
  arma::mat robust_matrix =
    solve(observed_fisher_info_matrix,
          trans(
            solve(observed_fisher_info_matrix,
                  meat_matrix)
          )
    );
  arma::vec lambda_vector = arma::zeros(params_no);
  for(int r = 1; r < params_no + 1; r++) {
    lambda_vector(r - 1) = -
      (
          trace(
            solve(observed_fisher_info_matrix,
                  partial_derivatives_matrix.rows((r - 1) * params_no,
                                                  r * params_no - 1)
            )
          ) +
            0.5 * trace(
                robust_matrix * second_derivatives_matrix.rows((r - 1) * params_no,
                                                               r * params_no - 1)
            )
      );
  }
  arma::vec ans =
    beta_vector + solve(observed_fisher_info_matrix, u_vector + lambda_vector);
  return ans;
}
//==============================================================================


//============================ update beta - jeffreys OR =======================
// [[Rcpp::export]]
arma::vec update_beta_jeffreys_or(const arma::vec & y_vector,
                                  const arma::mat & model_matrix,
                                  const arma::vec & id_vector,
                                  const arma::vec & repeated_vector,
                                  const arma::vec & weights_vector,
                                  const char * link,
                                  const arma::vec & beta_vector,
                                  const arma::vec & mu_vector,
                                  const arma::vec & eta_vector,
                                  const arma::vec & alpha_vector,
                                  const double & jeffreys_power){
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  int repeated_max = max(repeated_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::mat naive_matrix_inverse_derivative = arma::zeros(pow(params_no, 2), params_no);
  arma::mat meat_matrix = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec s_vector = y_vector - mu_vector;
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec mu_vector_i = mu_vector(id_vector_i);
    arma::vec weights_vector_i = weights_vector(id_vector_i);
    arma::mat d_matrix_i =
      arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max,
                                       alpha_vector);
    arma::mat v_matrix_inverse_i =
      get_v_matrix_inverse_or(mu_vector_i,
                              odds_ratios_vector_i,
                              weights_vector_i);
    arma::mat d_matrix_trans_v_matrix_inverse_i =
      trans(d_matrix_i) * v_matrix_inverse_i;
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector(id_vector_i);
    arma::mat v_matrix_inverse_delta_star_matrix_i =
      v_matrix_inverse_i * arma::diagmat(delta_star_vector(id_vector_i));
    naive_matrix_inverse_derivative +=
      trans(kron(d_matrix_i, d_matrix_i)) *
      (
          kronecker_sum_same(v_matrix_inverse_i *
            arma::diagmat(delta_star_vector(id_vector_i))) *
            kappa_matrix(mu_vector_i.n_elem) -
            kron(v_matrix_inverse_i, v_matrix_inverse_i) *
            get_v_matrix_mu_or(mu_vector_i,
                               odds_ratios_vector_i,
                               weights_vector_i)
      ) *
        d_matrix_i;
  }
  arma::vec lambda_vector =
    jeffreys_power *
    trans(naive_matrix_inverse_derivative) *
    vectorise(arma::inv(naive_matrix_inverse,
                        arma::inv_opts::allow_approx));
  arma::vec ans =
    beta_vector + solve(naive_matrix_inverse, u_vector + lambda_vector);
  return ans;
}


//============================ update beta - or ================================
// [[Rcpp::export]]
arma::vec update_beta_or(const arma::vec & y_vector,
                         const arma::mat & model_matrix,
                         const arma::vec & id_vector,
                         const arma::vec & repeated_vector,
                         const arma::vec & weights_vector,
                         const char * link,
                         const arma::vec & beta_vector,
                         const arma::vec & mu_vector,
                         const arma::vec & eta_vector,
                         const arma::vec & alpha_vector,
                         const double & jeffreys_power,
                         const char* method) {
  arma::vec ans;
  if(std::strcmp(method, "gee") == 0){
    ans = update_beta_gee_or(y_vector,
                             model_matrix,
                             id_vector,
                             repeated_vector,
                             weights_vector,
                             link,
                             beta_vector,
                             mu_vector,
                             eta_vector,
                             alpha_vector);
  }else if(std::strcmp(method, "brgee-naive") == 0){
    ans = update_beta_naive_or(y_vector,
                               model_matrix,
                               id_vector,
                               repeated_vector,
                               weights_vector,
                               link,
                               beta_vector,
                               mu_vector,
                               eta_vector,
                               alpha_vector);
  }else if(std::strcmp(method, "brgee-robust") == 0){
    ans = update_beta_robust_or(y_vector,
                                model_matrix,
                                id_vector,
                                repeated_vector,
                                weights_vector,
                                link,
                                beta_vector,
                                mu_vector,
                                eta_vector,
                                alpha_vector);
  }else if(std::strcmp(method, "brgee-empirical") == 0){
    ans = update_beta_empirical_or(y_vector,
                                   model_matrix,
                                   id_vector,
                                   repeated_vector,
                                   weights_vector,
                                   link,
                                   beta_vector,
                                   mu_vector,
                                   eta_vector,
                                   alpha_vector);
  }else if(std::strcmp(method, "pgee-jeffreys") == 0){
    ans = update_beta_jeffreys_or(y_vector,
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
  return(ans);
}
//==============================================================================


//=========================== fitting function =================================
// [[Rcpp::export]]
Rcpp::List fit_bingee_or(const arma::vec & y_vector,
                         const arma::mat & model_matrix,
                         const arma::vec & id_vector,
                         const arma::vec & repeated_vector,
                         const arma::vec & weight_vector,
                         const char * link,
                         arma::vec beta_vector,
                         const arma::vec & offset,
                         const int & maxiter,
                         const double & tolerance,
                         const int & step_maxiter,
                         const int & step_multiplier,
                         const double & jeffreys_power,
                         const char * method,
                         const arma::vec & alpha_vector) {
  int params_no = model_matrix.n_cols;
  arma::vec beta_vector_new = arma::zeros(params_no);
  arma::mat beta_hat_matrix = beta_vector;
  arma::vec stepsize_vector = arma::zeros(params_no);
  arma::vec criterion_vector = arma::zeros(maxiter);
  arma::vec beta_vector_inner = arma::zeros(params_no);
  arma::vec beta_vector_new_inner = arma::zeros(params_no);
  arma::vec stepsize_vector_inner = arma::zeros(params_no);
  arma::vec criterion_vector_inner = arma::zeros(maxiter);
  arma::vec eta_vector = model_matrix * beta_vector + offset;
  if(Rcpp::is_false(all(valideta(link, arma2vec(eta_vector)))))
    Rcpp::stop("invalid linear predictor - please another set of initial values for beta!!");
  arma::vec mu_vector = linkinv(link, arma2vec(eta_vector));
  if(Rcpp::is_false(all(validmu("binomial", arma2vec(mu_vector)))))
    Rcpp::stop("invalid fitted values - please another set of initial values for beta!!");
  for (int i = 1; i < maxiter + 1; i++){
    stepsize_vector = update_beta_or(y_vector,
                                     model_matrix,
                                     id_vector,
                                     repeated_vector,
                                     weight_vector,
                                     link,
                                     beta_vector,
                                     mu_vector,
                                     eta_vector,
                                     alpha_vector,
                                     jeffreys_power,
                                     method) - beta_vector;
    criterion_vector_inner(i - 1) = norm(stepsize_vector, "inf");
    beta_vector_inner = beta_vector;
    stepsize_vector_inner = stepsize_vector;
    for (int j = 1; j < step_maxiter + 1; j++){
      beta_vector_new_inner =
        beta_vector_inner + step_multiplier * pow(0.5, j - 1) * stepsize_vector_inner;
      eta_vector = model_matrix * beta_vector_new_inner + offset;
      if(Rcpp::is_false(all(valideta(link, arma2vec(eta_vector)))))
        Rcpp::stop("invalid initial linear predictor: please another set of initial values for beta!!");
      arma::vec mu_vector = linkinv(link, arma2vec(eta_vector));
      if(Rcpp::is_false(all(validmu("binomial", arma2vec(mu_vector)))))
        Rcpp::stop("invalid fitted values - please another set of initial values for beta!!");
      stepsize_vector_inner = update_beta_or(y_vector,
                                             model_matrix,
                                             id_vector,
                                             repeated_vector,
                                             weight_vector,
                                             link,
                                             beta_vector_new_inner,
                                             mu_vector,
                                             eta_vector,
                                             alpha_vector,
                                             jeffreys_power,
                                             method) - beta_vector_new_inner;
      beta_vector_new = beta_vector_new_inner;
      beta_vector_inner = beta_vector_new_inner;
      if (criterion_vector_inner(i - 1) > norm(stepsize_vector_inner, "inf")) break;
    }
    criterion_vector(i-1) = norm(stepsize_vector_inner, "inf");
    beta_vector = beta_vector_new;
    eta_vector = model_matrix * beta_vector_new + offset;
    if(Rcpp::is_false(all(valideta(link, arma2vec(eta_vector)))))
      Rcpp::stop("invalid linear predictor - please another set of initial values for beta!!");
    mu_vector = linkinv(link, arma2vec(eta_vector));
    if(Rcpp::is_false(all(validmu("binomial", arma2vec(mu_vector)))))
      Rcpp::stop("invalid fitted values - please another set of initial values for beta!!");
    beta_hat_matrix = join_rows(beta_hat_matrix, beta_vector_new);
    if(criterion_vector(i-1) <= tolerance) break;
  }
  Rcpp::List cov_matrices = get_covariance_matrices_or(y_vector,
                                                       model_matrix,
                                                       id_vector,
                                                       repeated_vector,
                                                       weight_vector,
                                                       link,
                                                       mu_vector,
                                                       eta_vector,
                                                       alpha_vector);
  Rcpp::List ans;
  ans["beta_hat"] = beta_vector;
  ans["beta_mat"] = beta_hat_matrix;
  ans["alpha"] = alpha_vector;
  ans["naive_covariance"] = cov_matrices[0];
  ans["robust_covariance"] = cov_matrices[1];
  ans["bc_covariance"] = cov_matrices[2];
  ans["criterion"] = criterion_vector;
  ans["criterion_ee"] = criterion_vector_inner;
  ans["eta"] = eta_vector;
  ans["residuals"] = y_vector - mu_vector;
  ans["fitted"] = mu_vector;
  return ans;
}
//==============================================================================
