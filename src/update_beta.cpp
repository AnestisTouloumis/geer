#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include "link_functions.h"
#include "variance_functions.h"
#include "nuisance_quantities.h"
#include "utils.h"


//============================ update beta - gee ===============================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_beta_gee(const arma::vec & y_vector,
                          const arma::mat & model_matrix,
                          const arma::vec & id_vector,
                          const arma::vec & repeated_vector,
                          const char* link,
                          const char* family,
                          const arma::vec & beta_vector,
                          const arma::vec & mu_vector,
                          const arma::vec & eta_vector,
                          const char * correlation_structure,
                          const arma::vec & alpha_vector,
                          const double & phi,
                          const arma::vec & weights_vector) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::mat correlation_matrix_inverse =
    get_correlation_matrix_inverse(correlation_structure,
                                   alpha_vector,
                                   max(repeated_vector));
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec s_vector = y_vector - mu_vector;
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::mat d_matrix_i = arma::diagmat(delta_vector(id_vector_i)) *
      model_matrix.rows(id_vector_i);
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(d_matrix_i) *
      get_weight_matrix_inverse(family,
                                mu_vector(id_vector_i),
                                repeated_vector(id_vector_i),
                                phi,
                                correlation_matrix_inverse,
                                weights_vector(id_vector_i));
    naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
    u_vector +=
      t_d_matrix_weight_matrix_inverse_i * s_vector(id_vector_i);
  }
  arma::vec ans = beta_vector + solve(naive_matrix_inverse, u_vector);
  return ans;
}
//==============================================================================


//============================ update beta - naive =============================
// [[Rcpp::export]]
arma::vec update_beta_naive(const arma::vec & y_vector,
                            const arma::mat & model_matrix,
                            const arma::vec & id_vector,
                            const arma::vec & repeated_vector,
                            const char* link,
                            const char* family,
                            const arma::vec & beta_vector,
                            const arma::vec & mu_vector,
                            const arma::vec & eta_vector,
                            const char * correlation_structure,
                            const arma::vec & alpha_vector,
                            const double & phi,
                            const arma::vec & weights_vector) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat gamma_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec alpha_star_vector =
    - 0.5 * variancemu(family, mu_vector) / variance(family, mu_vector);
    arma::vec s_vector = y_vector - mu_vector;
    arma::mat correlation_matrix_inverse =
      get_correlation_matrix_inverse(correlation_structure,
                                     alpha_vector,
                                     max(repeated_vector));
    for(int i=1; i < sample_size + 1; i++){
      arma::uvec id_vector_i = find(id_vector == i);
      arma::vec alpha_star_vector_i = alpha_star_vector(id_vector_i);
      arma::vec delta_star_vector_i = delta_star_vector(id_vector_i);
      arma::mat d_matrix_i =
        arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
      arma::mat weight_matrix_inverse_i =
        get_weight_matrix_inverse(family,
                                  mu_vector(id_vector_i),
                                  repeated_vector(id_vector_i),
                                  phi,
                                  correlation_matrix_inverse,
                                  weights_vector(id_vector_i));
      arma::mat t_d_matrix_weight_matrix_inverse_i =
        trans(d_matrix_i) * weight_matrix_inverse_i;
      naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
      u_vector += t_d_matrix_weight_matrix_inverse_i * s_vector(id_vector_i);
      arma::mat weight_matrix_inverse_alpha_plus_delta_star_matrix_i =
        weight_matrix_inverse_i *
        arma::diagmat(alpha_star_vector_i + delta_star_vector_i);
      gamma_matrix +=
        trans(kron(d_matrix_i, d_matrix_i)) *
        (kappa_right(trans(weight_matrix_inverse_alpha_plus_delta_star_matrix_i)) -
         kronecker_identity_right_kappa(
            weight_matrix_inverse_alpha_plus_delta_star_matrix_i
        ) -
          kronecker_left_identity_kappa(
            weight_matrix_inverse_i * arma::diagmat(delta_star_vector_i)
          )
        ) *
          d_matrix_i;
    }
    arma::vec gamma_vector = arma::zeros(params_no);
    for(int r = 1; r < params_no + 1; r++) {
      gamma_vector(r - 1) = - 0.5 *
        trace(solve(naive_matrix_inverse,
                    gamma_matrix.rows((r - 1) * params_no, r * params_no - 1)));
    }
    arma::vec ans = beta_vector + solve(naive_matrix_inverse, (u_vector + gamma_vector));
    return ans;
}
//==============================================================================


//============================ update beta - robust ============================
// [[Rcpp::export]]
arma::vec update_beta_robust(const arma::vec & y_vector,
                             const arma::mat & model_matrix,
                             const arma::vec & id_vector,
                             const arma::vec & repeated_vector,
                             const char* link,
                             const char* family,
                             const arma::vec & beta_vector,
                             const arma::vec & mu_vector,
                             const arma::vec & eta_vector,
                             const char * correlation_structure,
                             const arma::vec & alpha_vector,
                             const double & phi,
                             const arma::vec & weights_vector) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat gamma_matrix_one = arma::zeros(pow(params_no, 2), params_no);
  arma::mat gamma_matrix_two = arma::zeros(pow(params_no, 2), params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::mat meat_matrix = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2);
  arma::vec alpha_star_vector =
    - 0.5 * variancemu(family, mu_vector) / variance(family, mu_vector);
    arma::vec s_vector = y_vector - mu_vector;
    arma::mat correlation_matrix_inverse =
      get_correlation_matrix_inverse(correlation_structure,
                                     alpha_vector,
                                     max(repeated_vector));
    for(int i=1; i < sample_size + 1; i++){
      arma::uvec id_vector_i = find(id_vector == i);
      arma::vec s_vector_i = s_vector(id_vector_i);
      arma::vec alpha_star_vector_i = alpha_star_vector(id_vector_i);
      arma::mat d_matrix_i =
        arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
      arma::mat weight_matrix_inverse_i =
        get_weight_matrix_inverse(family,
                                  mu_vector(id_vector_i),
                                  repeated_vector(id_vector_i),
                                  phi,
                                  correlation_matrix_inverse,
                                  weights_vector(id_vector_i));
      arma::mat t_d_matrix_weight_matrix_inverse_i =
        trans(d_matrix_i) * weight_matrix_inverse_i;
      naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
      arma::vec u_vector_i = t_d_matrix_weight_matrix_inverse_i * s_vector_i;
      u_vector += u_vector_i;
      meat_matrix += u_vector_i * trans(u_vector_i);
      arma::mat weight_matrix_inverse_alpha_plus_delta_star_matrix_i =
        weight_matrix_inverse_i *
        arma::diagmat(alpha_star_vector_i + delta_star_vector(id_vector_i));
      arma::mat h_matrix_i1 =
        kappa_right(trans(weight_matrix_inverse_alpha_plus_delta_star_matrix_i));
      arma::mat weight_matrix_inverse_alpha_star_matrix_i =
        weight_matrix_inverse_i * arma::diagmat(alpha_star_vector_i);
      arma::mat kronecker_tdmatrix_tdmatrix = trans(kron(d_matrix_i, d_matrix_i));
      gamma_matrix_two -=
        kronecker_tdmatrix_tdmatrix *
        (h_matrix_i1 +
        kronecker_left_identity_kappa(
          weight_matrix_inverse_alpha_plus_delta_star_matrix_i +
            weight_matrix_inverse_alpha_star_matrix_i) +
            kronecker_identity_right_kappa(weight_matrix_inverse_alpha_plus_delta_star_matrix_i)) *
            d_matrix_i;
      gamma_matrix_one +=
        kronecker_tdmatrix_tdmatrix *
        (h_matrix_i1 +
        kronecker_left_identity_kappa(weight_matrix_inverse_alpha_star_matrix_i)) *
        s_vector_i * trans(u_vector_i);
    }
    arma::mat robust_matrix =
      solve(naive_matrix_inverse, trans(solve(naive_matrix_inverse, meat_matrix)));
    arma::vec gamma_vector = arma::zeros(params_no);
    for(int r = 1; r < params_no + 1; r++) {
      gamma_vector(r - 1) = -
        (trace(solve(naive_matrix_inverse,
                     gamma_matrix_one.rows((r - 1) * params_no, r * params_no - 1))) +
        0.5 * trace(robust_matrix * gamma_matrix_two.rows((r - 1) * params_no, r * params_no - 1))
        );
    }
    arma::vec ans = beta_vector + solve(naive_matrix_inverse, u_vector + gamma_vector);
    return ans;
}
//==============================================================================


//============================ update beta - empirical =========================
// [[Rcpp::export]]
arma::vec update_beta_empirical(const arma::vec & y_vector,
                                const arma::mat & model_matrix,
                                const arma::vec & id_vector,
                                const arma::vec & repeated_vector,
                                const char* link,
                                const char* family,
                                const arma::vec & beta_vector,
                                const arma::vec & mu_vector,
                                const arma::vec & eta_vector,
                                const char * correlation_structure,
                                const arma::vec & alpha_vector,
                                const double & phi,
                                const arma::vec & weights_vector){
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  arma::vec u_vector = arma::zeros(params_no);
  arma::mat gamma_matrix_one = arma::zeros(pow(params_no, 2), params_no);
  arma::mat gamma_matrix_two = arma::zeros(pow(params_no, 2), params_no);
  arma::mat gamma_matrix_three = arma::zeros(params_no, params_no);
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
  arma::mat correlation_matrix_inverse =
    get_correlation_matrix_inverse(correlation_structure,
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
    arma::mat weight_matrix_inverse_i =
      get_weight_matrix_inverse(family,
                                mu_vector(id_vector_i),
                                repeated_vector(id_vector_i),
                                phi,
                                correlation_matrix_inverse,
                                weights_vector(id_vector_i));
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(d_matrix_i) * weight_matrix_inverse_i;
    naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
    arma::vec u_vector_i = t_d_matrix_weight_matrix_inverse_i * s_vector_i;
    u_vector += u_vector_i;
    meat_matrix += u_vector_i * trans(u_vector_i);
    arma::vec w_vector_i = weight_matrix_inverse_i * s_vector_i;
    arma::mat epsilon_matrix_i =
      arma::diagmat(w_vector_i % alpha_star_plus_delta_star_vector_i) +
      weight_matrix_inverse_i *
      arma::diagmat(s_vector_i % alpha_star_vector_i - 1);
    arma::mat kronecker_tdmatrix_tdmatrix = trans(kron(d_matrix_i, d_matrix_i));
    gamma_matrix_one +=
      kronecker_tdmatrix_tdmatrix * arma::vectorise(trans(epsilon_matrix_i)) *
      trans(u_vector_i);
    gamma_matrix_two +=
      kronecker_tdmatrix_tdmatrix *
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
               weight_matrix_inverse_i *
               arma::diagmat(s_vector_i % alpha_tilde_star_vector_i -
               alpha_star_vector_i)
           )
      ) *
        d_matrix_i;
    gamma_matrix_three -= trans(d_matrix_i) * epsilon_matrix_i * d_matrix_i;
  }
  arma::mat gamma_three_matrix_inv = arma::inv(gamma_matrix_three, arma::inv_opts::allow_approx);
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


//============================ update beta - jeffreys ==========================
// [[Rcpp::export]]
arma::vec update_beta_jeffreys(const arma::vec & y_vector,
                               const arma::mat & model_matrix,
                               const arma::vec & id_vector,
                               const arma::vec & repeated_vector,
                               const char* link,
                               const char* family,
                               const arma::vec & beta_vector,
                               const arma::vec & mu_vector,
                               const arma::vec & eta_vector,
                               const char * correlation_structure,
                               const arma::vec & alpha_vector,
                               const double & phi,
                               const arma::vec & weights_vector,
                               const double & jeffreys_power) {
  int params_no = model_matrix.n_cols;
  int sample_size = max(id_vector);
  arma::mat u_vector = arma::zeros(params_no);
  arma::mat gamma_matrix = arma::zeros(pow(params_no, 2), params_no);
  arma::mat naive_matrix_inverse = arma::zeros(params_no, params_no);
  arma::vec delta_vector = mueta(link, eta_vector);
  arma::vec alpha_plus_delta_star_vector = mueta2(link, eta_vector)/pow(delta_vector, 2) -
    0.5 * variancemu(family, mu_vector) / variance(family, mu_vector);
  arma::vec s_vector = y_vector - mu_vector;
  arma::mat correlation_matrix_inverse =
    get_correlation_matrix_inverse(correlation_structure,
                                   alpha_vector,
                                   max(repeated_vector));
  for(int i=1; i < sample_size + 1; i++){
    arma::uvec id_vector_i = find(id_vector == i);
    arma::mat d_matrix_i =
      arma::diagmat(delta_vector(id_vector_i)) * model_matrix.rows(id_vector_i);
    arma::mat weight_matrix_inverse_i =
      get_weight_matrix_inverse(family,
                                mu_vector(id_vector_i),
                                repeated_vector(id_vector_i),
                                phi,
                                correlation_matrix_inverse,
                                weights_vector(id_vector_i));
    arma::mat t_d_matrix_weight_matrix_inverse_i =
      trans(d_matrix_i) * weight_matrix_inverse_i;
    naive_matrix_inverse += t_d_matrix_weight_matrix_inverse_i * d_matrix_i;
    u_vector += t_d_matrix_weight_matrix_inverse_i * s_vector(id_vector_i);
    arma::mat weight_matrix_inverse_alpha_plus_delta_star_matrix_i =
      weight_matrix_inverse_i * arma::diagmat(alpha_plus_delta_star_vector(id_vector_i));
    gamma_matrix +=
      trans(kron(d_matrix_i, d_matrix_i)) *
      (kronecker_left_identity_kappa(weight_matrix_inverse_alpha_plus_delta_star_matrix_i) +
      kronecker_identity_right_kappa(weight_matrix_inverse_alpha_plus_delta_star_matrix_i)) *
      d_matrix_i;
  }
  //  arma::vec gamma_vector = arma::zeros(params_no);
  //  for(int r = 1; r < params_no + 1; r++) {
  //    gamma_vector(r - 1) = 0.5 *
  //    trace(solve(naive_matrix_inverse,
  //                gamma_matrix.rows((r - 1) * params_no, r * params_no - 1)));
  //}
  arma::vec gamma_vector = jeffreys_power * trans(gamma_matrix) * vectorise(inv(naive_matrix_inverse));
  arma::vec ans = beta_vector + solve(naive_matrix_inverse, u_vector + gamma_vector);
  return ans;
}
//==============================================================================




//============================ update beta =====================================
// [[Rcpp::export]]
arma::vec update_beta(const arma::vec & y_vector,
                      const arma::mat & model_matrix,
                      const arma::vec & id_vector,
                      const arma::vec & repeated_vector,
                      const char* link,
                      const char* family,
                      const arma::vec & beta_vector,
                      const arma::vec & mu_vector,
                      const arma::vec & eta_vector,
                      const char * correlation_structure,
                      const arma::vec & alpha_vector,
                      const double & phi,
                      const char* type,
                      const arma::vec & weights_vector,
                      const double & jeffreys_power) {
  arma::vec ans(beta_vector.n_elem);
  if(std::strcmp(type, "gee") == 0){
    ans = update_beta_gee(y_vector,
                          model_matrix,
                          id_vector,
                          repeated_vector,
                          link,
                          family,
                          beta_vector,
                          mu_vector,
                          eta_vector,
                          correlation_structure,
                          alpha_vector,
                          phi,
                          weights_vector);
  }else if(std::strcmp(type, "brgee_naive") == 0){
    ans = update_beta_naive(y_vector,
                            model_matrix,
                            id_vector,
                            repeated_vector,
                            link,
                            family,
                            beta_vector,
                            mu_vector,
                            eta_vector,
                            correlation_structure,
                            alpha_vector,
                            phi,
                            weights_vector);
  }else if(std::strcmp(type, "brgee_robust") == 0){
    ans = update_beta_robust(y_vector,
                             model_matrix,
                             id_vector,
                             repeated_vector,
                             link,
                             family,
                             beta_vector,
                             mu_vector,
                             eta_vector,
                             correlation_structure,
                             alpha_vector,
                             phi,
                             weights_vector);
  }else if(std::strcmp(type, "brgee_empirical") == 0){
    ans = update_beta_empirical(y_vector,
                                model_matrix,
                                id_vector,
                                repeated_vector,
                                link,
                                family,
                                beta_vector,
                                mu_vector,
                                eta_vector,
                                correlation_structure,
                                alpha_vector,
                                phi,
                                weights_vector);
  }else if(std::strcmp(type, "pgee_jeffreys") == 0){
    ans = update_beta_jeffreys(y_vector,
                               model_matrix,
                               id_vector,
                               repeated_vector,
                               link,
                               family,
                               beta_vector,
                               mu_vector,
                               eta_vector,
                               correlation_structure,
                               alpha_vector,
                               phi,
                               weights_vector,
                               jeffreys_power);
  }
  return(ans);
}
//==============================================================================
