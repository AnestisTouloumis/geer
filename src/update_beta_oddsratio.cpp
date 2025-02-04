#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link_functions.h"
#include "utils.h"
#include "nuisance_quantities_oddsratio.h"


//============================ update beta - gee OR ============================
// [[Rcpp::export]]
arma::vec update_beta_gee_or(const arma::vec & y_vector,
                              const arma::mat & model_matrix,
                              const arma::vec & id_vector,
                              const arma::vec & repeated_vector,
                              const char * link,
                              const arma::vec & beta_vector,
                              const arma::vec & mu_vector,
                              const arma::vec & eta_vector,
                              const arma::vec & alpha_vector,
                              const arma::vec & weights_vector) {
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
      arma::diagmat(delta_vector(id_vector_i)) *
      model_matrix.rows(id_vector_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max,
                                       alpha_vector);
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(d_matrix_i) *
      get_weight_matrix_inverse_or(mu_vector(id_vector_i),
                                   odds_ratios_vector_i,
                                   weights_vector(id_vector_i));
    naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
    u_vector +=
      t_d_matrix_weight_matrix_inverse_i * s_vector(id_vector_i);
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
                               const char * link,
                               const arma::vec & beta_vector,
                               const arma::vec & mu_vector,
                               const arma::vec & eta_vector,
                               const arma::vec & alpha_vector,
                               const arma::vec & weights_vector) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  int repeated_max = max(repeated_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat gamma_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec s_vector = y_vector - mu_vector;
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    double cluster_size_i = id_vector_i.n_elem;
    arma::vec mu_vector_i = mu_vector(id_vector_i);
    arma::vec weights_vector_i = weights_vector(id_vector_i);
    arma::mat d_matrix_i =
      diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max,
                                       alpha_vector);
    arma::mat weight_matrix_inverse_i =
      get_weight_matrix_inverse_or(mu_vector_i, odds_ratios_vector_i, weights_vector_i);
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(d_matrix_i) * weight_matrix_inverse_i;
    naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
    u_vector +=
      t_d_matrix_weight_matrix_inverse_i * s_vector(id_vector_i);
    arma::mat delta_star_matrix_i = arma::diagmat(delta_star_vector(id_vector_i));
    arma::mat g_matrix_i = get_g_matrix(mu_vector_i, odds_ratios_vector_i);
    arma::mat identity_matrix_i = arma::eye(cluster_size_i, cluster_size_i);
    arma::mat weight_matrix_inverse_delta_star_matrix_i =
      weight_matrix_inverse_i * delta_star_matrix_i;
    arma::mat weight_matrix_inverse_trans_g_matrix_i =
      weight_matrix_inverse_i * trans(g_matrix_i);
    arma::mat kappa_matrix_delta_matrix_i = kappa_right(delta_star_matrix_i);
    arma::mat vectorised_weight_matrix_inverse_trans_mu_vector_i =
      arma::vectorise(weight_matrix_inverse_i) * trans(mu_vector_i);
    arma::vec weight_matrix_inverse_mu_vector_i =
      weight_matrix_inverse_i * mu_vector_i;
    arma::mat gamma_matrix_two =
      (kronecker_left_identity_kappa(weight_matrix_inverse_delta_star_matrix_i) +
      kronecker_identity_right_kappa(weight_matrix_inverse_delta_star_matrix_i)) -
      (kronecker_identity_right_kappa(weight_matrix_inverse_i) * g_matrix_i +
      kronecker_left_identity_kappa(weight_matrix_inverse_trans_g_matrix_i +
      weight_matrix_inverse_i) -
      kappa_matrix_delta_matrix_i -
      vectorised_weight_matrix_inverse_trans_mu_vector_i) * weight_matrix_inverse_i -
      (kron(weight_matrix_inverse_i, weight_matrix_inverse_i)) *
      get_weight_matrix_mu_or(mu_vector_i, odds_ratios_vector_i, weights_vector_i) +
      kronecker_vector_matrix(weight_matrix_inverse_mu_vector_i, weight_matrix_inverse_i);
    arma::mat gamma_matrix_one =
      (kappa_matrix_delta_matrix_i -
      kronecker_left_identity_kappa(weight_matrix_inverse_i) *
      (g_matrix_i + arma::eye(cluster_size_i, cluster_size_i)) -
      kronecker_left_identity_kappa(weight_matrix_inverse_trans_g_matrix_i) +
      vectorised_weight_matrix_inverse_trans_mu_vector_i +
      kronecker_vector_identity(weight_matrix_inverse_mu_vector_i)) *
      weight_matrix_inverse_i;
    gamma_matrix -= trans(kron(d_matrix_i, d_matrix_i)) *
      (gamma_matrix_one - 0.5 * gamma_matrix_two) *
      d_matrix_i;
  }
  arma::vec gamma_vector = arma::zeros(params_no);
  for(int r = 1; r < params_no + 1; r++) {
    gamma_vector(r - 1) =
      trace(solve(naive_matrix_inverse,
                  gamma_matrix.rows((r - 1) * params_no, r * params_no - 1)));
  }
  arma::vec ans = beta_vector +
    solve(naive_matrix_inverse, u_vector + gamma_vector);
  return ans;
}
//==============================================================================


//============================ update beta - robust OR =========================
// [[Rcpp::export]]
arma::vec update_beta_robust_or(const arma::vec & y_vector,
                                const arma::mat & model_matrix,
                                const arma::vec & id_vector,
                                const arma::vec & repeated_vector,
                                const char * link,
                                const arma::vec & beta_vector,
                                const arma::vec & mu_vector,
                                const arma::vec & eta_vector,
                                const arma::vec & alpha_vector,
                                const arma::vec & weights_vector) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  int repeated_max = max(repeated_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat gamma_matrix_one = arma::zeros(pow(params_no, 2), params_no);
  arma::mat gamma_matrix_two = arma::zeros(pow(params_no, 2), params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::mat meat_matrix = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec s_vector = y_vector - mu_vector;
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    double cluster_size_i = id_vector_i.n_elem;
    arma::vec mu_vector_i = mu_vector(id_vector_i);
    arma::vec weights_vector_i = weights_vector(id_vector_i);
    arma::mat d_matrix_i =
      diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max,
                                       alpha_vector);
    arma::mat weight_matrix_inverse_i =
      get_weight_matrix_inverse_or(mu_vector_i, odds_ratios_vector_i, weights_vector_i);
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(d_matrix_i) * weight_matrix_inverse_i;
    naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
    arma::vec s_vector_i = s_vector(id_vector_i);
    arma::vec u_vector_i = t_d_matrix_weight_matrix_inverse_i * s_vector_i;
    u_vector += u_vector_i;
    meat_matrix += u_vector_i * trans(u_vector_i);
    arma::mat delta_star_matrix_i =
      arma::diagmat(delta_star_vector(id_vector_i));
    arma::mat g_matrix_i = get_g_matrix(mu_vector_i, odds_ratios_vector_i);
    arma::mat weight_matrix_inverse_delta_star_matrix_i =
      weight_matrix_inverse_i * delta_star_matrix_i;
    arma::mat vectorised_weight_matrix_inverse_trans_mu_vector_i =
      arma::vectorise(weight_matrix_inverse_i) * trans(mu_vector_i);
    arma::mat weight_matrix_inverse_trans_g_matrix_i =
      weight_matrix_inverse_i * trans(g_matrix_i);
    arma::mat kappa_matrix_delta_matrix_i = kappa_right(delta_star_matrix_i);
    arma::vec weight_matrix_inverse_mu_vector_i =
      weight_matrix_inverse_i * mu_vector_i;
    arma::mat kronecker_tdmatrix_tdmatrix = trans(kron(d_matrix_i, d_matrix_i));
    gamma_matrix_two -=
      kronecker_tdmatrix_tdmatrix * (
          (kronecker_left_identity_kappa(weight_matrix_inverse_delta_star_matrix_i) +
            kronecker_identity_right_kappa(weight_matrix_inverse_delta_star_matrix_i)) -
            (kronecker_identity_right_kappa(weight_matrix_inverse_i) * g_matrix_i +
            kronecker_left_identity_kappa(weight_matrix_inverse_trans_g_matrix_i +
            weight_matrix_inverse_i) -
            kappa_matrix_delta_matrix_i -
            vectorised_weight_matrix_inverse_trans_mu_vector_i) * weight_matrix_inverse_i -
            (kron(weight_matrix_inverse_i, weight_matrix_inverse_i)) *
            get_weight_matrix_mu_or(mu_vector_i, odds_ratios_vector_i, weights_vector_i) +
            kronecker_vector_matrix(weight_matrix_inverse_mu_vector_i, weight_matrix_inverse_i)) *
            d_matrix_i;
    gamma_matrix_one +=
      kronecker_tdmatrix_tdmatrix *
      (kappa_matrix_delta_matrix_i -
      kronecker_left_identity_kappa(weight_matrix_inverse_i) *
      (g_matrix_i + arma::eye(cluster_size_i, cluster_size_i)) -
      kronecker_left_identity_kappa(weight_matrix_inverse_trans_g_matrix_i) +
      vectorised_weight_matrix_inverse_trans_mu_vector_i +
      kronecker_vector_identity(weight_matrix_inverse_mu_vector_i)) *
      weight_matrix_inverse_i *
      s_vector_i *
      trans(u_vector_i);
  }
  arma::mat robust_matrix =
    solve(naive_matrix_inverse, trans(solve(naive_matrix_inverse, meat_matrix)));
  arma::vec gamma_vector = arma::zeros(params_no);
  for(int r = 1; r < params_no + 1; r++) {
    gamma_vector(r - 1) = -
      (trace(solve(naive_matrix_inverse, gamma_matrix_one.rows((r - 1) * params_no, r * params_no - 1))) +
      0.5 * trace(robust_matrix * gamma_matrix_two.rows((r - 1) * params_no, r * params_no - 1))
      );
  }
  arma::vec ans = beta_vector + solve(naive_matrix_inverse, u_vector + gamma_vector);
  return ans;
}
//==============================================================================


//============================ update beta - empirical OR ======================
// [[Rcpp::export]]
arma::vec update_beta_empirical_or(const arma::vec & y_vector,
                                   const arma::mat & model_matrix,
                                   const arma::vec & id_vector,
                                   const arma::vec & repeated_vector,
                                   const char * link,
                                   const arma::vec & beta_vector,
                                   const arma::vec & mu_vector,
                                   const arma::vec & eta_vector,
                                   const arma::vec & alpha_vector,
                                   const arma::vec & weights_vector) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  int repeated_max = max(repeated_vector);
  arma::vec u_vector = arma::zeros(params_no);
  arma::mat gamma_matrix_one = arma::zeros(pow(params_no, 2), params_no);
  arma::mat gamma_matrix_two = arma::zeros(pow(params_no, 2), params_no);
  arma::mat gamma_matrix_three = arma::zeros(params_no, params_no);
  arma::mat meat_matrix = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec s_vector = y_vector - mu_vector;
  arma::vec delta_tilde_star_vector =
    (delta_vector % mueta3(link, eta_vector) - 2 * pow(mueta2(link, eta_vector), 2))/
      pow(delta_vector, 4);
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    double cluster_size_i = id_vector_i.n_elem;
    arma::vec mu_vector_i = mu_vector(id_vector_i);
    arma::vec eta_vector_i = eta_vector(id_vector_i);
    arma::vec mueta_vector_i = mueta(link, eta_vector_i);
    arma::vec weights_vector_i = weights_vector(id_vector_i);
    arma::mat d_matrix_i =
      arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max,
                                       alpha_vector);
    arma::mat weight_matrix_inverse_i =
      get_weight_matrix_inverse_or(mu_vector_i, odds_ratios_vector_i, weights_vector_i);
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(d_matrix_i) * weight_matrix_inverse_i;
    arma::vec s_vector_i = s_vector(id_vector_i);
    arma::vec u_vector_i = t_d_matrix_weight_matrix_inverse_i * s_vector_i;
    u_vector += u_vector_i;
    meat_matrix += u_vector_i * trans(u_vector_i);
    arma::vec weight_matrix_inverse_s_vector_i =
      weight_matrix_inverse_i * s_vector_i;
    arma::mat w_matrix_i = arma::diagmat(weight_matrix_inverse_s_vector_i);
    arma::mat delta_star_matrix_i = arma::diagmat(delta_star_vector(id_vector_i));
    arma::mat g_matrix_i = get_g_matrix(mu_vector_i, odds_ratios_vector_i);
    arma::mat l_matrix_i =
      arma::diagmat(g_matrix_i * weight_matrix_inverse_s_vector_i);
    arma::mat identity_matrix_i = arma::eye(cluster_size_i, cluster_size_i);
    arma::vec weight_matrix_inverse_mu_vector_i = weight_matrix_inverse_i * mu_vector_i;
    double help_vec = as_scalar(trans(s_vector_i) * weight_matrix_inverse_mu_vector_i);
    arma::mat epsilon_matrix_i =
      (delta_star_matrix_i -
      weight_matrix_inverse_i * (identity_matrix_i + trans(g_matrix_i))) *
      w_matrix_i +
      help_vec * weight_matrix_inverse_i -
      weight_matrix_inverse_i *
      (l_matrix_i +
      identity_matrix_i -
      mu_vector_i * trans(s_vector_i) * weight_matrix_inverse_i);
    // first row
    arma::mat transpose_epsilon_matrix_deriv_i1 =
      kappa_right(w_matrix_i * arma::diagmat(delta_tilde_star_vector(id_vector_i))) +
      arma::vectorise(weight_matrix_inverse_i) * trans(weight_matrix_inverse_s_vector_i) -
      kronecker_vector_matrix(weight_matrix_inverse_mu_vector_i, weight_matrix_inverse_i) +
      kron(weight_matrix_inverse_i, weight_matrix_inverse_s_vector_i);
    //second row
    arma::mat transpose_epsilon_matrix_deriv_i2 =
      kappa_right(delta_star_matrix_i) +
      arma::vectorise(weight_matrix_inverse_i) * trans(mu_vector_i) -
      kronecker_left_identity_kappa(weight_matrix_inverse_i) * g_matrix_i -
      kronecker_left_identity_kappa(weight_matrix_inverse_i * (trans(g_matrix_i) + identity_matrix_i));
    transpose_epsilon_matrix_deriv_i2 =
      transpose_epsilon_matrix_deriv_i2 * (epsilon_matrix_i - w_matrix_i * delta_star_matrix_i);
    // derivative of v_matrix inverse
    arma::mat weight_matrix_inverse_deriv_i =
      - kron(weight_matrix_inverse_i, weight_matrix_inverse_i) *
      get_weight_matrix_mu_or(mu_vector_i, odds_ratios_vector_i, weights_vector_i);
    // third row
    arma::mat transpose_epsilon_matrix_deriv_i3 =
      kron(weight_matrix_inverse_mu_vector_i * trans(s_vector_i) - identity_matrix_i,
           identity_matrix_i) +
             kron(identity_matrix_i,
                  weight_matrix_inverse_s_vector_i * trans(mu_vector_i) -
                    w_matrix_i * (identity_matrix_i + g_matrix_i) - l_matrix_i) +
                    help_vec * kron(identity_matrix_i, identity_matrix_i);
    transpose_epsilon_matrix_deriv_i3 =
      transpose_epsilon_matrix_deriv_i3 * weight_matrix_inverse_deriv_i;
    // fourth row
    arma::mat transpose_epsilon_matrix_deriv_i4 =
      (kron(weight_matrix_inverse_i, w_matrix_i) +
      kronecker_left_identity_kappa(weight_matrix_inverse_i) *
      trans(kronecker_vector_identity(weight_matrix_inverse_s_vector_i))) *
      get_g_matrix_mu(mu_vector_i, odds_ratios_vector_i);
    arma::mat transpose_epsilon_matrix_deriv_i =
      transpose_epsilon_matrix_deriv_i1 +
      transpose_epsilon_matrix_deriv_i2 +
      transpose_epsilon_matrix_deriv_i3 -
      transpose_epsilon_matrix_deriv_i4;
    arma::mat h_matrix_i1 =
      kronecker_left_identity_kappa(epsilon_matrix_i * delta_star_matrix_i) +
      kronecker_identity_right_kappa(trans(epsilon_matrix_i) * delta_star_matrix_i);
    arma::mat tdmatrix_epsilon_dmatrix = trans(d_matrix_i) * epsilon_matrix_i * d_matrix_i;
    gamma_matrix_one += arma::vectorise(trans(tdmatrix_epsilon_dmatrix)) * trans(u_vector_i);
    gamma_matrix_two +=
      trans(kron(d_matrix_i, d_matrix_i)) * (h_matrix_i1 + transpose_epsilon_matrix_deriv_i) * d_matrix_i;
    gamma_matrix_three -= tdmatrix_epsilon_dmatrix;
  }
  arma::mat gamma_three_matrix_inv = arma::pinv(gamma_matrix_three);
  arma::mat robust_matrix = gamma_three_matrix_inv * meat_matrix *
    trans(gamma_three_matrix_inv);
  arma::vec gamma_vector = arma::zeros(params_no);
  for(int r = 1; r < params_no + 1; r++) {
    gamma_vector(r - 1) = -
      (trace(gamma_three_matrix_inv * gamma_matrix_one.rows((r - 1) * params_no, r * params_no - 1)) +
      0.5 * trace(robust_matrix * gamma_matrix_two.rows((r - 1) * params_no, r * params_no - 1))
      );
  }
  arma::vec ans = beta_vector + gamma_three_matrix_inv * (u_vector + gamma_vector);
  return ans;
  }
//==============================================================================


//============================ update beta - jeffreys OR =======================
// [[Rcpp::export]]
arma::vec update_beta_jeffreys_or(const arma::vec & y_vector,
                                  const arma::mat & model_matrix,
                                  const arma::vec & id_vector,
                                  const arma::vec & repeated_vector,
                                  const char * link,
                                  const arma::vec & beta_vector,
                                  const arma::vec & mu_vector,
                                  const arma::vec & eta_vector,
                                  const arma::vec & alpha_vector,
                                  const arma::vec & weights_vector,
                                  const double & jeffreys_power){
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  int repeated_max = max(repeated_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat gamma_matrix = arma::zeros(pow(params_no, 2), params_no);
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
      arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::vec odds_ratios_vector_i =
      get_subject_specific_odds_ratios(repeated_vector(id_vector_i),
                                       repeated_max, alpha_vector);
    arma::mat weight_matrix_inverse_i =
      get_weight_matrix_inverse_or(mu_vector_i, odds_ratios_vector_i, weights_vector_i);
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(d_matrix_i) * weight_matrix_inverse_i;
    naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
    u_vector +=
      t_d_matrix_weight_matrix_inverse_i * s_vector(id_vector_i);
    arma::mat weight_matrix_inverse_delta_star_matrix_i =
      weight_matrix_inverse_i * arma::diagmat(delta_star_vector(id_vector_i));
    gamma_matrix +=
      trans(kron(d_matrix_i, d_matrix_i)) *
      (kronecker_sum_same(weight_matrix_inverse_delta_star_matrix_i) * kappa_matrix(mu_vector_i.n_elem) -
      kron(weight_matrix_inverse_i, weight_matrix_inverse_i) *
      get_weight_matrix_mu_or(mu_vector_i, odds_ratios_vector_i, weights_vector_i)) *
      d_matrix_i;
  }
//  arma::vec gamma_vector = arma::zeros(params_no);
//  for(int r = 1; r < params_no + 1; r++) {
//    gamma_vector(r - 1) = 0.5 *
//      trace(solve(naive_matrix_inverse,
//                  gamma_matrix.rows((r - 1) * params_no, r * params_no - 1)));
//  }
  arma::vec gamma_vector = jeffreys_power * trans(gamma_matrix) * vectorise(pinv(naive_matrix_inverse));
  arma::vec ans = beta_vector + solve(naive_matrix_inverse, u_vector + gamma_vector);
  return ans;
}


//==============================================================================
//============================ update beta - or ================================
// [[Rcpp::export]]
arma::vec update_beta_or(const arma::vec & y_vector,
                         const arma::mat & model_matrix,
                         const arma::vec & id_vector,
                         const arma::vec & repeated_vector,
                         const char * link,
                         const arma::vec & beta_vector,
                         const arma::vec & mu_vector,
                         const arma::vec & eta_vector,
                         const arma::vec & alpha_vector,
                         const char* type,
                         const arma::vec & weights_vector,
                         const double & jeffreys_power) {
  arma::vec ans;
  if(std::strcmp(type, "gee") == 0){
    ans = update_beta_gee_or(y_vector,
                             model_matrix,
                             id_vector,
                             repeated_vector,
                             link,
                             beta_vector,
                             mu_vector,
                             eta_vector,
                             alpha_vector,
                             weights_vector);
  }else if(std::strcmp(type, "brgee_naive") == 0){
    ans = update_beta_naive_or(y_vector,
                               model_matrix,
                               id_vector,
                               repeated_vector,
                               link,
                               beta_vector,
                               mu_vector,
                               eta_vector,
                               alpha_vector,
                               weights_vector);
  }else if(std::strcmp(type, "brgee_robust") == 0){
    ans = update_beta_robust_or(y_vector,
                                model_matrix,
                                id_vector,
                                repeated_vector,
                                link,
                                beta_vector,
                                mu_vector,
                                eta_vector,
                                alpha_vector,
                                weights_vector);
  }else if(std::strcmp(type, "brgee_empirical") == 0){
    ans = update_beta_empirical_or(y_vector,
                                   model_matrix,
                                   id_vector,
                                   repeated_vector,
                                   link,
                                   beta_vector,
                                   mu_vector,
                                   eta_vector,
                                   alpha_vector,
                                   weights_vector);
  }else if(std::strcmp(type, "pgee_jeffreys") == 0){
    ans = update_beta_jeffreys_or(y_vector,
                                  model_matrix,
                                  id_vector,
                                  repeated_vector,
                                  link,
                                  beta_vector,
                                  mu_vector,
                                  eta_vector,
                                  alpha_vector,
                                  weights_vector,
                                  jeffreys_power);
  }
  return(ans);
}
//==============================================================================
