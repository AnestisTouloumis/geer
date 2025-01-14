#include <RcppArmadillo.h>
#include "nuisance_quantities.h"
#include "link_functions.h"
#include "covariance_matrices.h"
#include "variance_functions.h"
#include "update_beta.h"


//=========================== fitting function =================================
// [[Rcpp::export]]
Rcpp::List fit_geesolver(const arma::vec & y_vector,
                         const arma::mat & model_matrix,
                         const arma::vec & id_vector,
                         const arma::vec & repeated_vector,
                         const char * link,
                         const char * family,
                         arma::vec beta_vector,
                         const char * correlation_structure,
                         const int & mdependence,
                         int use_params,
                         const int & maxiter,
                         const double & tolerance,
                         const arma::vec & offset,
                         double phi,
                         const int & phi_fixed,
                         arma::vec alpha_vector,
                         const int & alpha_fixed,
                         const char * type,
                         const arma::vec & weights_vector,
                         const int & step_maxiter,
                         const int & step_multiplier,
                         const double & jeffreys_power) {
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
    Rcpp::Rcerr << "invalid linear predictor\n";
  arma::vec mu_vector = linkinv(link, arma2vec(eta_vector));
  if(Rcpp::is_false(all(Rcpp::is_finite(arma2vec(mu_vector)) &
     (arma2vec(mu_vector) > 0) & (arma2vec(mu_vector) < 1))))
    Rcpp::Rcerr << "invalid fitted values\n";
  arma::vec pearson_residuals_vector =
    get_pearson_residuals(family,
                          y_vector,
                          mu_vector);
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
    stepsize_vector = update_beta(y_vector,
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
                                  type,
                                  weights_vector,
                                  jeffreys_power) - beta_vector;
    criterion_vector_inner(i - 1) = norm(stepsize_vector, "inf");
    beta_vector_inner = beta_vector;
    stepsize_vector_inner = stepsize_vector;
    for (int j = 1; j < step_maxiter + 1; j++){
      beta_vector_new_inner =
        beta_vector_inner + step_multiplier * pow(0.5, j - 1) * stepsize_vector_inner;
      eta_vector = model_matrix * beta_vector_new_inner + offset;
      if(Rcpp::is_false(all(valideta(link, arma2vec(eta_vector)))))
        Rcpp::Rcerr << "invalid linear predictor\n";
      mu_vector = linkinv(link, arma2vec(eta_vector));
      if(Rcpp::is_false(all(Rcpp::is_finite(arma2vec(mu_vector)) &
         (arma2vec(mu_vector) > 0) & (arma2vec(mu_vector) < 1))))
        Rcpp::Rcerr << "invalid fitted values\n";
      pearson_residuals_vector = get_pearson_residuals(family,
                                                       y_vector,
                                                       mu_vector);
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
      stepsize_vector_inner = update_beta(y_vector,
                                          model_matrix,
                                          id_vector,
                                          repeated_vector,
                                          link,
                                          family,
                                          beta_vector_new_inner,
                                          mu_vector,
                                          eta_vector,
                                          correlation_structure,
                                          alpha_vector,
                                          phi,
                                          type,
                                          weights_vector,
                                          jeffreys_power) - beta_vector_new_inner;
      beta_vector_new = beta_vector_new_inner;
      beta_vector_inner = beta_vector_new_inner;
      if (criterion_vector_inner(i - 1) > norm(stepsize_vector_inner, "inf")) break;
    }
    criterion_vector(i-1) = norm(stepsize_vector_inner, "inf");
    beta_vector = beta_vector_new;
    eta_vector = model_matrix * beta_vector_new + offset;
    if(Rcpp::is_false(all(valideta(link, arma2vec(eta_vector)))))
      Rcpp::Rcerr << "invalid linear predictor\n";
    mu_vector = linkinv(link, arma2vec(eta_vector));
    if(Rcpp::is_false(all(Rcpp::is_finite(arma2vec(mu_vector)) &
       (arma2vec(mu_vector) > 0) & (arma2vec(mu_vector) < 1))))
      Rcpp::Rcerr << "invalid fitted values\n";
    beta_hat_matrix = join_rows(beta_hat_matrix, beta_vector_new);
    pearson_residuals_vector = get_pearson_residuals(family,
                                                     y_vector,
                                                     mu_vector);
    if(phi_fixed == 0)
      phi = get_phi_hat(pearson_residuals_vector, use_params);
    if(criterion_vector(i-1) <= tolerance) break;
  }
  Rcpp::List cov_matrices = get_covariance_matrices(y_vector,
                                                    model_matrix,
                                                    id_vector,
                                                    repeated_vector,
                                                    link,
                                                    family,
                                                    mu_vector,
                                                    eta_vector,
                                                    correlation_structure,
                                                    alpha_vector,
                                                    phi,
                                                    weights_vector);
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
  return ans;
}
//==============================================================================
