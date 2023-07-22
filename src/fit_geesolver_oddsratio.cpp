#include <RcppArmadillo.h>
#include "link_functions.h"
#include "covariance_matrices_oddsratio.h"
#include "update_beta_oddsratio.h"


//=========================== fitting function =================================
// [[Rcpp::export]]
Rcpp::List fit_bingee_or(const arma::vec & y_vector,
                         const arma::mat & model_matrix,
                         const arma::vec & id_vector,
                         const arma::vec & repeated_vector,
                         const char * link,
                         arma::vec beta_vector,
                         const int & maxiter,
                         const double & tolerance,
                         const arma::vec & offset,
                         const arma::vec & alpha_vector,
                         const char * type) {
  int params_no = model_matrix.n_cols;
  arma::vec criterion_vector = arma::zeros(maxiter);
  arma::vec beta_vector_new = arma::zeros(params_no);
  arma::mat beta_hat_matrix = beta_vector;
  arma::vec eta_vector = model_matrix * beta_vector + offset;
  if(Rcpp::is_false(all(valideta(link, arma2vec(eta_vector)))))
    Rcpp::Rcerr << "invalid linear predictor\n";
  arma::vec mu_vector = linkinv(link, arma2vec(eta_vector));
  if(Rcpp::is_false(all(Rcpp::is_finite(arma2vec(mu_vector)) &
     (arma2vec(mu_vector) > 0) & (arma2vec(mu_vector) < 1))))
    Rcpp::Rcerr << "invalid fitted values\n";
  for (int i = 1; i < maxiter + 1; i++){
    beta_vector_new = update_beta_or(y_vector,
                                     model_matrix,
                                     id_vector,
                                     repeated_vector,
                                     link,
                                     beta_vector,
                                     mu_vector,
                                     eta_vector,
                                     alpha_vector,
                                     type);
    criterion_vector(i-1) =
      max(abs(beta_vector - beta_vector_new));
    beta_vector = beta_vector_new;
    beta_hat_matrix = join_rows(beta_hat_matrix, beta_vector_new);
    eta_vector = model_matrix * beta_vector + offset;
    if(Rcpp::is_false(all(valideta(link, arma2vec(eta_vector)))))
      Rcpp::Rcerr << "invalid linear predictor\n";
    mu_vector = linkinv(link, arma2vec(eta_vector));
    if(Rcpp::is_false(all(Rcpp::is_finite(arma2vec(mu_vector)) &
       (arma2vec(mu_vector) > 0) & (arma2vec(mu_vector) < 1))))
      Rcpp::Rcerr << "invalid fitted values\n";
    if(criterion_vector(i-1) <= tolerance) break;
  }
  Rcpp::List cov_matrices = get_covariance_matrices_or(y_vector,
                                                       model_matrix,
                                                       id_vector,
                                                       repeated_vector,
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
  ans["eta"] = eta_vector;
  ans["residuals"] = y_vector - mu_vector;
  ans["fitted"] = mu_vector;
  return ans;
}
//==============================================================================
