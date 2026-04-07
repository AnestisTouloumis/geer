#include <RcppArmadillo.h>
using namespace Rcpp;

//============================ arma to vec =====================================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector arma2vec(const arma::vec & x) {
  Rcpp::NumericVector ans = Rcpp::NumericVector(x.begin(), x.end());
  return(ans);
  }

//==============================================================================


//============================ vec to arma =====================================
// [[Rcpp::export]]
arma::vec vec2arma(const Rcpp::NumericVector & x) {
  arma::vec ans = Rcpp::as<arma::vec>(x);
  return(ans);
  }
//==============================================================================

//============================ subset matrix x[1:y, 1:y] =======================
// [[Rcpp::depends(RcppArmadillo)]]
arma::mat subset_matrix(const arma::mat& x, const arma::vec& y) {
 // arma::uvec z = arma::conv_to<arma::uvec>::from(arma::vec(y));
//  return x.submat(z - 1, z - 1);
  arma::uvec z = arma::conv_to<arma::uvec>::from(y);
  z -= 1;
  return x.submat(z, z);
}
//==============================================================================


//============================ kappa matrix ====================================
arma::mat kappa_matrix(int dimension) {
  //arma::mat ans = arma::zeros(dimension * dimension, dimension);
  //for(int j = 1; j < dimension + 1; j++) {
  //  ans((j - 1) * dimension + j - 1, j - 1) = 1;
  //}
  //return ans;
  arma::mat ans(dimension * dimension, dimension, arma::fill::zeros);
  for (arma::uword j = 0; j < dimension; ++j) {
    ans(j * dimension + j, j) = 1.0;
  }
  return ans;
}

//==============================================================================


//============================ kronecker(X, I) * kappa =========================
arma::mat kronecker_left_identity_kappa(const arma::mat& x) {
  //const int dimension = x.n_rows;
  //const int dim_sq = dimension * dimension;
  //arma::mat ans = arma::zeros(dim_sq, dimension);
  //for(int j = 1; j < dimension + 1; j++) {
  //  ans.rows((j - 1) * dimension, j * dimension - 1).diag() =
  //    x.row(j - 1);
  //}
  //return ans;
  const arma::uword dimension = x.n_rows;
  arma::mat ans(dimension * dimension, dimension, arma::fill::zeros);
  for (arma::uword j = 0; j < dimension; ++j) {
    const arma::uword base = j * dimension;
    for (arma::uword i = 0; i < dimension; ++i) {
      ans(base + i, i) = x(j, i);
    }
  }
  return ans;
}
//==============================================================================


//============================ kronecker(I, X) * kappa =========================
arma::mat kronecker_identity_right_kappa(const arma::mat& x) {
 // const int dimension = x.n_rows;
  //const int dim_sq = dimension * dimension;
  //arma::mat ans = arma::zeros(dim_sq, dimension);
  //for(int  j = 1; j < dimension + 1; j++) {
  //  ans.submat((j - 1) * dimension, j - 1, j * dimension - 1, j - 1) =
  //    x.col(j - 1);
  //}
  //return(ans);
  const arma::uword dimension = x.n_rows;
  arma::mat ans(dimension * dimension, dimension, arma::fill::zeros);
  for (arma::uword j = 0; j < dimension; ++j) {
    ans.submat(j * dimension, j, (j + 1) * dimension - 1, j) = x.col(j);
  }
  return ans;
}
//==============================================================================


//============================ Kappa * X  ======================================
arma::mat kappa_right(const arma::mat& x) {
//  const int dimension = x.n_rows;
//  const int dim_sq = dimension * dimension;
//  arma::mat ans = arma::zeros(dim_sq, dimension);
//  for(int  j = 1; j < dimension + 1; j++) {
//    ans.row((j - 1) * (dimension + 1)) = x.row(j - 1);
//  }
//  return(ans);
const arma::uword dimension = x.n_rows;
  arma::mat ans(dimension * dimension, dimension, arma::fill::zeros);
  arma::uvec rows = arma::regspace<arma::uvec>(0, dimension - 1) * (dimension + 1);
  ans.rows(rows) = x;
  return ans;
}
//==============================================================================


//============================ kronecker direct sum of the same matrix =========
arma::mat kronecker_sum_same(const arma::mat& x) {
  const arma::uword dimension = x.n_rows;
  const arma::mat identity_mat = arma::eye(dimension, dimension);
  return arma::kron(x, identity_mat) + arma::kron(identity_mat, x);
}
//==============================================================================


//============================ kronecker(x, I)  ================================
arma::mat kronecker_vector_identity(const arma::vec& x) {
  const arma::uword dimension = x.n_rows;
  return arma::kron(x, arma::eye(dimension, dimension));
}
//==============================================================================


//============================ kronecker(x, Y)   ===============================
arma::mat kronecker_vector_matrix(const arma::vec& x, const arma::mat& y) {
  return arma::kron(x, y);
}
//==============================================================================


//============================ X^{-1} y  =======================================
arma::vec solve_chol_or_lu_vec(const arma::mat& X, const arma::vec& y) {
  arma::mat chol_upper;
  if (arma::chol(chol_upper, X)) {
    arma::vec chol_intermediate = arma::solve(arma::trimatl(chol_upper.t()), y);
    return arma::solve(arma::trimatu(chol_upper), chol_intermediate);
  }
  arma::mat lu_lower, lu_upper, permutation_matrix;
  arma::lu(lu_lower, lu_upper, permutation_matrix, X);
  arma::vec permutation_vector = permutation_matrix * y;
  arma::vec lu_forward = arma::solve(arma::trimatl(lu_lower), permutation_vector, arma::solve_opts::allow_ugly);
  return arma::solve(arma::trimatu(lu_upper), lu_forward,  arma::solve_opts::allow_ugly);
}

//=========================== X^{-1} Y  ========================================
arma::mat solve_chol_or_lu_mat(const arma::mat& X, const arma::mat& Y) {
  arma::mat chol_upper;
  if (arma::chol(chol_upper, X)) {
    arma::mat Z = arma::solve(arma::trimatl(chol_upper.t()), Y);
    return arma::solve(arma::trimatu(chol_upper), Z);
  }
  arma::mat lu_lower, lu_upper, permutation_matrix;
  arma::lu(lu_lower, lu_upper, permutation_matrix, X);
  arma::mat permuted_matrix = permutation_matrix * Y;
  arma::mat Z = arma::solve(arma::trimatl(lu_lower), permuted_matrix, arma::solve_opts::allow_ugly);
  return arma::solve(arma::trimatu(lu_upper), Z,  arma::solve_opts::allow_ugly);
}

// ---- lambda blocks
arma::vec lambda_from_blocks_chol_or_lu(const arma::mat& A,
                                        const arma::mat& lambda_matrix) {
  const arma::uword p = A.n_rows;
  arma::vec lambda_vector(p, arma::fill::zeros);

  arma::mat R;
  if (arma::chol(R, A)) {
    for (arma::uword r = 0; r < p; ++r) {
      const arma::mat block = lambda_matrix.rows(r * p, (r + 1) * p - 1);
      arma::mat Y = arma::solve(arma::trimatl(R.t()), block);
      arma::mat X = arma::solve(arma::trimatu(R), Y);
      lambda_vector[r] = -0.5 * arma::trace(X);
    }
    return lambda_vector;
  }

  arma::mat L, U, P;
  arma::lu(L, U, P, A);

  for (arma::uword r = 0; r < p; ++r) {
    const arma::mat block = lambda_matrix.rows(r * p, (r + 1) * p - 1);
    const arma::mat PB = P * block;
    const arma::mat Y  = arma::solve(arma::trimatl(L), PB, arma::solve_opts::allow_ugly);
    const arma::mat X  = arma::solve(arma::trimatu(U), Y,  arma::solve_opts::allow_ugly);
    lambda_vector[r] = -0.5 * arma::trace(X);
  }
  return lambda_vector;
}
