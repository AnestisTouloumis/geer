#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "nuisance_quantities_cc.h"
#include "utils.h"
#include "link_functions.h"
#include "variance_functions.h"
#include "covariance_matrices.h"

//============================ update beta - gee ===============================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_beta_gee_cc(const arma::vec & y_vector,
                             const arma::mat & model_matrix,
                             const arma::vec & id_vector,
                             const arma::vec & repeated_vector,
                             const arma::vec & weights_vector,
                             const char* link,
                             const char* family,
                             const arma::vec & beta_vector,
                             const arma::vec & mu_vector,
                             const arma::vec & eta_vector,
                             const char * correlation_structure,
                             const arma::vec & alpha_vector,
                             const double & phi) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  int repeated_max = max(repeated_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure,
                           alpha_vector,
                           repeated_max);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec s_vector = y_vector - mu_vector;
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::mat d_matrix_i =
      arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::mat d_matrix_trans_v_matrix_inverse_i =
      trans(
        solve(
          get_v_matrix_cc(family,
                          mu_vector(id_vector_i),
                          repeated_vector(id_vector_i),
                          phi,
                          correlation_matrix,
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


//============================ update beta - naive =============================
// [[Rcpp::export]]
arma::vec update_beta_naive_cc(const arma::vec & y_vector,
                               const arma::mat & model_matrix,
                               const arma::vec & id_vector,
                               const arma::vec & repeated_vector,
                               const arma::vec & weights_vector,
                               const char* link,
                               const char* family,
                               const arma::vec & beta_vector,
                               const arma::vec & mu_vector,
                               const arma::vec & eta_vector,
                               const char * correlation_structure,
                               const arma::vec & alpha_vector,
                               const double & phi) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat lambda_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec alpha_star_vector =
    - 0.5 * variancemu(family, mu_vector) / variance(family, mu_vector);
    arma::vec s_vector = y_vector - mu_vector;
    arma::mat correlation_matrix =
      get_correlation_matrix(correlation_structure,
                             alpha_vector,
                             max(repeated_vector));
    for(int i=1; i < sample_size + 1; i++){
      arma::uvec id_vector_i = find(id_vector == i);
      arma::vec alpha_star_vector_i = alpha_star_vector(id_vector_i);
      arma::vec delta_star_vector_i = delta_star_vector(id_vector_i);
      arma::mat d_matrix_i =
        arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
      arma::mat v_matrix_inverse_i =
        arma::inv(get_v_matrix_cc(family,
                                  mu_vector(id_vector_i),
                                  repeated_vector(id_vector_i),
                                  phi,
                                  correlation_matrix,
                                  weights_vector(id_vector_i)),
                                  arma::inv_opts::allow_approx);
      arma::mat d_matrix_trans_v_matrix_inverse_i =
        trans(d_matrix_i) * v_matrix_inverse_i;
      naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector(id_vector_i);
      arma::mat v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i =
        v_matrix_inverse_i *
        arma::diagmat(alpha_star_vector_i + delta_star_vector_i);
      lambda_matrix +=
        trans(kron(d_matrix_i, d_matrix_i)) *
        (kappa_right(trans(v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i)) -
        kronecker_identity_right_kappa(
          v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i
        ) -
          kronecker_left_identity_kappa(
            v_matrix_inverse_i * arma::diagmat(delta_star_vector_i)
          )
        ) *
          d_matrix_i;
    }
    arma::vec lambda_vector = arma::zeros(params_no);
    for(int r = 1; r < params_no + 1; r++) {
      lambda_vector(r - 1) = - 0.5 *
        trace(solve(naive_matrix_inverse,
                    lambda_matrix.rows((r - 1) * params_no, r * params_no - 1)));
    }
    arma::vec ans =
      beta_vector + solve(naive_matrix_inverse, (u_vector + lambda_vector));
    return ans;
}
//==============================================================================


//============================ update beta - robust ============================
// [[Rcpp::export]]
arma::vec update_beta_robust_cc(const arma::vec & y_vector,
                                const arma::mat & model_matrix,
                                const arma::vec & id_vector,
                                const arma::vec & repeated_vector,
                                const arma::vec & weights_vector,
                                const char* link,
                                const char* family,
                                const arma::vec & beta_vector,
                                const arma::vec & mu_vector,
                                const arma::vec & eta_vector,
                                const char * correlation_structure,
                                const arma::vec & alpha_vector,
                                const double & phi) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat partial_derivatives_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat second_derivatives_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::mat meat_matrix = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec alpha_star_vector =
    - 0.5 * variancemu(family, mu_vector) / variance(family, mu_vector);
    arma::vec s_vector = y_vector - mu_vector;
    arma::mat correlation_matrix =
      get_correlation_matrix(correlation_structure,
                             alpha_vector,
                             max(repeated_vector));
    for(int i=1; i < sample_size + 1; i++){
      arma::uvec id_vector_i = find(id_vector == i);
      arma::vec s_vector_i = s_vector(id_vector_i);
      arma::vec alpha_star_vector_i = alpha_star_vector(id_vector_i);
      arma::mat d_matrix_i =
        arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
      arma::mat v_matrix_inverse_i =
        arma::inv(get_v_matrix_cc(family,
                                  mu_vector(id_vector_i),
                                  repeated_vector(id_vector_i),
                                  phi,
                                  correlation_matrix,
                                  weights_vector(id_vector_i)),
                                  arma::inv_opts::allow_approx);
      arma::mat d_matrix_trans_v_matrix_inverse_i =
        trans(d_matrix_i) * v_matrix_inverse_i;
      naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      arma::vec u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector_i;
      u_vector += u_vector_i;
      meat_matrix += u_vector_i * trans(u_vector_i);
      arma::mat v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i =
        v_matrix_inverse_i *
        arma::diagmat(alpha_star_vector_i + delta_star_vector(id_vector_i));
      arma::mat kappa_matrix_delta_star_matrix_plus_alpha_star_matrix_v_matrix_inverse_i =
        kappa_right(trans(v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i));
      arma::mat v_matrix_inverse_alpha_star_matrix_i =
        v_matrix_inverse_i * arma::diagmat(alpha_star_vector_i);
      arma::mat kron_d_matrix_trans_d_matrix_trans_i = trans(kron(d_matrix_i, d_matrix_i));
      second_derivatives_matrix -=
        kron_d_matrix_trans_d_matrix_trans_i *
        (kappa_matrix_delta_star_matrix_plus_alpha_star_matrix_v_matrix_inverse_i +
        kronecker_left_identity_kappa(
          v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i +
            v_matrix_inverse_alpha_star_matrix_i) +
            kronecker_identity_right_kappa(
              v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i
            )) *
              d_matrix_i;
      partial_derivatives_matrix +=
        kron_d_matrix_trans_d_matrix_trans_i *
        (kappa_matrix_delta_star_matrix_plus_alpha_star_matrix_v_matrix_inverse_i +
        kronecker_left_identity_kappa(v_matrix_inverse_alpha_star_matrix_i)) *
        s_vector_i * trans(u_vector_i);
    }
    arma::mat robust_matrix =
      solve(naive_matrix_inverse, trans(solve(naive_matrix_inverse, meat_matrix)));
    arma::vec lambda_vector = arma::zeros(params_no);
    for(int r = 1; r < params_no + 1; r++) {
      lambda_vector(r - 1) = -
        (trace(
            solve(naive_matrix_inverse,
                  partial_derivatives_matrix.rows((r - 1) * params_no,
                                                  r * params_no - 1))) +
                                                    0.5 * trace(
                                                        robust_matrix * second_derivatives_matrix.rows((r - 1) * params_no,
                                                                                                       r * params_no - 1))
        );
    }
    arma::vec ans = beta_vector + solve(naive_matrix_inverse, u_vector + lambda_vector);
    return ans;
}
//==============================================================================


//============================ update beta - empirical =========================
// [[Rcpp::export]]
arma::vec update_beta_empirical_cc(const arma::vec & y_vector,
                                   const arma::mat & model_matrix,
                                   const arma::vec & id_vector,
                                   const arma::vec & repeated_vector,
                                   const arma::vec & weights_vector,
                                   const char* link,
                                   const char* family,
                                   const arma::vec & beta_vector,
                                   const arma::vec & mu_vector,
                                   const arma::vec & eta_vector,
                                   const char * correlation_structure,
                                   const arma::vec & alpha_vector,
                                   const double & phi){
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  arma::vec u_vector = arma::zeros(params_no);
  arma::mat partial_derivatives_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat second_derivatives_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat observed_fisher_info_matrix = arma::zeros(params_no, params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::mat meat_matrix = arma::zeros(params_no, params_no);
  arma::vec s_vector = y_vector - mu_vector;
  arma::vec mueta2_vector = mueta2(link, eta_vector);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2_vector/pow(delta_vector, 2);
  arma::vec delta_tilde_star_vector =
    (delta_vector % mueta3(link, eta_vector) - 2 * pow(mueta2_vector, 2))/
      pow(delta_vector, 4);
  arma::vec variance_vector = variance(family, mu_vector);
  arma::vec variancemu_vector = variancemu(family, mu_vector);
  arma::vec alpha_star_vector = - 0.5 * variancemu_vector / variance_vector;
  arma::vec alpha_tilde_star_vector =
    0.5 * (pow(variancemu_vector, 2)/variance_vector - variancemu2(family, mu_vector))/
      variance_vector;
  arma::vec alpha_star_plus_delta_star_vector = alpha_star_vector + delta_star_vector;
  arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure,
                           alpha_vector,
                           max(repeated_vector));
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec s_vector_i = s_vector(id_vector_i);
    arma::vec alpha_star_vector_i = alpha_star_vector(id_vector_i);
    arma::vec alpha_tilde_star_vector_i = alpha_tilde_star_vector(id_vector_i);
    arma::vec alpha_star_plus_delta_star_vector_i =
      alpha_star_plus_delta_star_vector(id_vector_i);
    arma::mat d_matrix_i =
      diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::mat v_matrix_inverse_i =
      arma::inv(get_v_matrix_cc(family,
                                mu_vector(id_vector_i),
                                repeated_vector(id_vector_i),
                                phi,
                                correlation_matrix,
                                weights_vector(id_vector_i)),
                                arma::inv_opts::allow_approx);
    arma::mat d_matrix_trans_v_matrix_inverse_i =
      trans(d_matrix_i) * v_matrix_inverse_i;
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    arma::vec u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector_i;
    u_vector += u_vector_i;
    meat_matrix += u_vector_i * trans(u_vector_i);
    arma::vec w_vector_i = v_matrix_inverse_i * s_vector_i;
    arma::mat epsilon_matrix_i =
      arma::diagmat(w_vector_i % alpha_star_plus_delta_star_vector_i) +
      v_matrix_inverse_i *
      arma::diagmat(s_vector_i % alpha_star_vector_i - 1);
    arma::mat observed_fisher_info_matrix_i = trans(d_matrix_i) * epsilon_matrix_i * d_matrix_i;
    observed_fisher_info_matrix -= observed_fisher_info_matrix_i;
    partial_derivatives_matrix +=
      arma::vectorise(trans(observed_fisher_info_matrix_i)) * trans(u_vector_i);
    second_derivatives_matrix +=
      trans(kron(d_matrix_i, d_matrix_i)) *
      (kappa_right(
          arma::diagmat(alpha_star_plus_delta_star_vector_i) * epsilon_matrix_i -
            arma::diagmat(alpha_star_plus_delta_star_vector_i %
            (alpha_star_plus_delta_star_vector_i + alpha_star_vector_i) %
            w_vector_i) +
            arma::diagmat((alpha_tilde_star_vector_i + delta_tilde_star_vector(id_vector_i)) %
            w_vector_i))
         + kronecker_identity_right_kappa(
             trans(epsilon_matrix_i) * arma::diagmat(alpha_star_plus_delta_star_vector_i)
         ) +
           kronecker_left_identity_kappa(
             epsilon_matrix_i *
               arma::diagmat(alpha_star_plus_delta_star_vector_i) +
               v_matrix_inverse_i *
               arma::diagmat(s_vector_i % alpha_tilde_star_vector_i -
               alpha_star_vector_i)
           )
      ) *
        d_matrix_i;
  }
  arma::mat robust_matrix =
    solve(observed_fisher_info_matrix,
          trans(solve(observed_fisher_info_matrix, meat_matrix)));
  arma::vec lambda_vector = arma::zeros(params_no);
  for(int r = 1; r < params_no + 1; r++) {
    lambda_vector(r - 1) = -
      (trace(solve(observed_fisher_info_matrix,
                   partial_derivatives_matrix.rows((r - 1) * params_no,
                                                   r * params_no - 1))) +
                                                     0.5 * trace(robust_matrix * second_derivatives_matrix.rows((r - 1) * params_no,
                                                                                                                r * params_no - 1))
      );
  }
  arma::vec ans =
    beta_vector + solve(observed_fisher_info_matrix, u_vector + lambda_vector);
  return ans;
}
//==============================================================================


//============================ update beta - jeffreys ==========================
// [[Rcpp::export]]
arma::vec update_beta_jeffreys_cc(const arma::vec & y_vector,
                                  const arma::mat & model_matrix,
                                  const arma::vec & id_vector,
                                  const arma::vec & repeated_vector,
                                  const arma::vec & weights_vector,
                                  const char* link,
                                  const char* family,
                                  const arma::vec & beta_vector,
                                  const arma::vec & mu_vector,
                                  const arma::vec & eta_vector,
                                  const char * correlation_structure,
                                  const arma::vec & alpha_vector,
                                  const double & phi,
                                  const double & jeffreys_power) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat t_naive_matrix_inverse_derivative = arma::zeros(pow(params_no, 2), params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec alpha_plus_delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2) -
    0.5 * variancemu(family, mu_vector) / variance(family, mu_vector);
  arma::vec s_vector = y_vector - mu_vector;
  arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure,
                           alpha_vector,
                           max(repeated_vector));
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::mat d_matrix_i =
      arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::mat v_matrix_inverse_i =
      arma::inv(get_v_matrix_cc(family,
                                mu_vector(id_vector_i),
                                repeated_vector(id_vector_i),
                                phi,
                                correlation_matrix,
                                weights_vector(id_vector_i)),
                                arma::inv_opts::allow_approx);
    arma::mat d_matrix_trans_v_matrix_inverse_i =
      trans(d_matrix_i) * v_matrix_inverse_i;
    naive_matrix_inverse += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
    u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector(id_vector_i);
    arma::mat v_matrix_inverse_alpha_matrix_plus_delta_star_matrix_i =
      v_matrix_inverse_i * arma::diagmat(alpha_plus_delta_star_vector(id_vector_i));
    t_naive_matrix_inverse_derivative +=
      trans(kron(d_matrix_i, d_matrix_i)) *
      (kronecker_left_identity_kappa(v_matrix_inverse_alpha_matrix_plus_delta_star_matrix_i) +
      kronecker_identity_right_kappa(v_matrix_inverse_alpha_matrix_plus_delta_star_matrix_i)) *
      d_matrix_i;
  }
  arma::vec lambda_vector =
    jeffreys_power *
    trans(t_naive_matrix_inverse_derivative) *
    vectorise(arma::inv(naive_matrix_inverse,
                        arma::inv_opts::allow_approx));
  arma::vec ans = beta_vector + solve(naive_matrix_inverse, u_vector + lambda_vector);
  return ans;
}
//==============================================================================



//============================ update beta =====================================
// [[Rcpp::export]]
arma::vec update_beta_cc(const arma::vec & y_vector,
                         const arma::mat & model_matrix,
                         const arma::vec & id_vector,
                         const arma::vec & repeated_vector,
                         const arma::vec & weights_vector,
                         const char* link,
                         const char* family,
                         const arma::vec & beta_vector,
                         const arma::vec & mu_vector,
                         const arma::vec & eta_vector,
                         const char * correlation_structure,
                         const arma::vec & alpha_vector,
                         const double & phi,
                         const double & jeffreys_power,
                         const char* method) {
  arma::vec ans(beta_vector.n_elem);
  if(std::strcmp(method, "gee") == 0){
    ans = update_beta_gee_cc(y_vector,
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
  }else if(std::strcmp(method, "brgee_naive") == 0){
    ans = update_beta_naive_cc(y_vector,
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
  }else if(std::strcmp(method, "brgee_robust") == 0){
    ans = update_beta_robust_cc(y_vector,
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
  }else if(std::strcmp(method, "brgee_empirical") == 0){
    ans = update_beta_empirical_cc(y_vector,
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
  }else if(std::strcmp(method, "pgee_jeffreys") == 0){
    ans = update_beta_jeffreys_cc(y_vector,
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
  return(ans);
}
//==============================================================================

//=========================== fitting function =================================
// [[Rcpp::export]]
Rcpp::List fit_geesolver_cc(const arma::vec & y_vector,
                            const arma::mat & model_matrix,
                            const arma::vec & id_vector,
                            const arma::vec & repeated_vector,
                            const arma::vec & weights_vector,
                            const char * link,
                            const char * family,
                            arma::vec beta_vector,
                            const arma::vec & offset,
                            const int & maxiter,
                            const double & tolerance,
                            const int & step_maxiter,
                            const int & step_multiplier,
                            const double & jeffreys_power,
                            const char * method,
                            int use_params,
                            arma::vec alpha_vector,
                            const int & alpha_fixed,
                            const char * correlation_structure,
                            const int & mdependence,
                            double phi,
                            const int & phi_fixed) {
  int params_no = model_matrix.n_cols;
  use_params = use_params * params_no;
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
    Rcpp::stop("invalid initial linear predictor: please another set of initial values for beta!!");
  arma::vec mu_vector = linkinv(link, arma2vec(eta_vector));
  if(Rcpp::is_false(all(validmu(family, arma2vec(mu_vector)))))
    Rcpp::stop("invalid initial fitted values: please another set of initial values for beta!!");
  arma::vec pearson_residuals_vector =
    get_pearson_residuals(family,
                          y_vector,
                          mu_vector,
                          weights_vector);
  if(phi_fixed == 0)
    phi = get_phi_hat(pearson_residuals_vector, use_params);
  for (int i=1; i<maxiter+1; i++){
    if(alpha_fixed == 0)
      alpha_vector = get_alpha_hat(correlation_structure,
                                   pearson_residuals_vector,
                                   id_vector,
                                   repeated_vector,
                                   phi,
                                   use_params,
                                   mdependence);
    stepsize_vector = update_beta_cc(y_vector,
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
    criterion_vector_inner(i - 1) = norm(stepsize_vector, "inf");
    beta_vector_inner = beta_vector;
    stepsize_vector_inner = stepsize_vector;
    for (int j = 1; j < step_maxiter + 1; j++){
      beta_vector_new_inner =
        beta_vector_inner + step_multiplier * pow(0.5, j - 1) * stepsize_vector_inner;
      eta_vector = model_matrix * beta_vector_new_inner + offset;
      if(Rcpp::is_false(all(valideta(link, arma2vec(eta_vector)))))
        Rcpp::stop("invalid linear predictor - please another set of initial values for beta!!");
      mu_vector = linkinv(link, arma2vec(eta_vector));
      if(Rcpp::is_false(all(validmu(family, arma2vec(mu_vector)))))
        Rcpp::stop("invalid fitted values - please another set of initial values for beta!!");
      pearson_residuals_vector = get_pearson_residuals(family,
                                                       y_vector,
                                                       mu_vector,
                                                       weights_vector);
      if(phi_fixed == 0)
        phi = get_phi_hat(pearson_residuals_vector, use_params);
      if(alpha_fixed == 0)
        alpha_vector = get_alpha_hat(correlation_structure,
                                     pearson_residuals_vector,
                                     id_vector,
                                     repeated_vector,
                                     phi,
                                     use_params,
                                     mdependence);
      stepsize_vector_inner = update_beta_cc(y_vector,
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
      if (criterion_vector_inner(i - 1) > norm(stepsize_vector_inner, "inf")) break;
    }
    criterion_vector(i-1) = norm(stepsize_vector_inner, "inf");
    beta_vector = beta_vector_new;
    eta_vector = model_matrix * beta_vector_new + offset;
    if(Rcpp::is_false(all(valideta(link, arma2vec(eta_vector)))))
      Rcpp::stop("invalid linear predictor - please another set of initial values for beta!!");
    mu_vector = linkinv(link, arma2vec(eta_vector));
    if(Rcpp::is_false(all(validmu(family, arma2vec(mu_vector)))))
      Rcpp::stop("invalid fitted values - please another set of initial values for beta!!");
    beta_hat_matrix = join_rows(beta_hat_matrix, beta_vector_new);
    pearson_residuals_vector = get_pearson_residuals(family,
                                                     y_vector,
                                                     mu_vector,
                                                     weights_vector);
    if(phi_fixed == 0)
      phi = get_phi_hat(pearson_residuals_vector, use_params);
    if(criterion_vector(i-1) <= tolerance) break;
  }
  Rcpp::List cov_matrices = get_covariance_matrices_cc(y_vector,
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
