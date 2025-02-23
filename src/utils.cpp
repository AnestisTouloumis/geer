#include <RcppArmadillo.h>
using namespace Rcpp;


//============================ subset matrix x[1:y, 1:y] =======================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat subset_matrix(arma::mat x, arma::vec y) {
  arma::uvec z = arma::conv_to<arma::uvec>::from(arma::vec(y));
  return(x.submat(z - 1, z - 1));
}
//==============================================================================


//============================ kappa matrix ====================================
// [[Rcpp::export]]
arma::mat kappa_matrix(int dimension) {
  arma::mat ans = arma::zeros(pow(dimension, 2), dimension);
  for(int j = 1; j < dimension + 1; j++) {
    ans((j - 1) * dimension + j - 1, j - 1) = 1;
  }
  return(ans);
}
//==============================================================================


//============================ kronecker direct sum of the same matrix =========
// [[Rcpp::export]]
arma::mat kronecker_sum_same(arma::mat x) {
  arma::mat identity_matrix = arma::eye(x.n_rows, x.n_rows);
  arma::mat ans = kron(x, identity_matrix) + kron(identity_matrix, x);
  return(ans);
}
//==============================================================================


//============================ kronecker(X, I) * kappa =========================
// [[Rcpp::export]]
arma::mat kronecker_left_identity_kappa(arma::mat x) {
  int dimension = x.n_rows;
  arma::mat ans = arma::zeros(pow(dimension, 2), dimension);
  for(int j = 1; j < dimension + 1; j++) {
    ans.rows((j - 1) * dimension, j * dimension - 1).diag() =
      x.row(j - 1);
  }
  return(ans);
}
//==============================================================================


//============================ kronecker(I, X) * kappa =========================
// [[Rcpp::export]]
arma::mat kronecker_identity_right_kappa(arma::mat x) {
  int dimension = x.n_rows;
  arma::mat ans = arma::zeros(pow(dimension, 2), dimension);
  for(int  j = 1; j < dimension + 1; j++) {
    ans.submat((j - 1) * dimension, j - 1, j * dimension - 1, j - 1) =
      x.col(j - 1);
  }
  return(ans);
}
//==============================================================================


//============================ Kappa * X  ======================================
// [[Rcpp::export]]
arma::mat kappa_right(arma::mat x) {
  int dimension = x.n_rows;
  arma::mat ans = arma::zeros(pow(dimension, 2), dimension);
  for(int  j = 1; j < dimension + 1; j++) {
    ans.row((j - 1) * (dimension + 1)) = x.row(j - 1);
  }
  return(ans);
}
//==============================================================================


//============================ kronecker(x, I)  ================================
// [[Rcpp::export]]
arma::mat kronecker_vector_identity(arma::vec x) {
  int dimension = x.n_rows;
  arma::mat ans = arma::zeros(pow(dimension, 2), dimension);
  for(int  j = 1; j < dimension + 1; j++) {
    ans.rows((j - 1) * dimension, j * dimension - 1).diag() += x(j - 1);
  }
  return(ans);
}
//==============================================================================


//============================ kronecker(x, Y)   ===============================
// [[Rcpp::export]]
arma::mat kronecker_vector_matrix(arma::vec x, arma::mat y) {
  int dimension = x.n_rows;
  arma::mat ans = arma::repmat(y, dimension, 1);
  for(int  j = 1; j < dimension + 1; j++) {
    ans.rows((j - 1) * dimension, j * dimension - 1) *= x(j - 1);
  }
  return(ans);
}
//==============================================================================


//============================ kronecker    ===============================
// [[Rcpp::export]]
arma::mat kronecker_rcpp(arma::mat x, arma::mat y) {
  arma::mat ans = kron(x, y);
  return(ans);
}
//==============================================================================
