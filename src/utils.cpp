#include "utils.h"
#include <cmath>


//============================ arma to vec =====================================
Rcpp::NumericVector arma2vec(const arma::vec& x) {
  return Rcpp::NumericVector(x.begin(), x.end());
}
//==============================================================================


//============================ vec to arma =====================================
arma::vec vec2arma(const Rcpp::NumericVector& x) {
  return Rcpp::as<arma::vec>(x);
}
//==============================================================================


//============================ subset matrix x[y, y] ===========================
arma::mat subset_matrix(const arma::mat& x, const arma::vec& y) {
  arma::uvec z = arma::conv_to<arma::uvec>::from(y);
  z -= 1;
  return x.submat(z, z);
}
//==============================================================================


//============================ kappa matrix ====================================
arma::mat kappa_matrix(const arma::uword dimension) {
  arma::mat ans(dimension * dimension, dimension, arma::fill::zeros);

  for (arma::uword j = 0; j < dimension; ++j) {
    ans(j * dimension + j, j) = 1.0;
  }
  return ans;
}
//==============================================================================


//============================ kronecker(X, I) * kappa =========================
arma::mat kronecker_left_identity_kappa(const arma::mat& x) {
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
  const arma::uword dimension = x.n_rows;
  arma::mat ans(dimension * dimension, dimension, arma::fill::zeros);

  for (arma::uword j = 0; j < dimension; ++j) {
    ans.submat(j * dimension, j, (j + 1) * dimension - 1, j) = x.col(j);
  }

  return ans;
}
//==============================================================================


//============================ Kappa * X =======================================
arma::mat kappa_right(const arma::mat& x) {
  const arma::uword dimension = x.n_rows;
  arma::mat ans(dimension * dimension, dimension, arma::fill::zeros);

  const arma::uvec rows =
    arma::regspace<arma::uvec>(0, dimension - 1) * (dimension + 1);
  ans.rows(rows) = x;

  return ans;
}
//==============================================================================


//============================ kronecker direct sum of same matrix =============
arma::mat kronecker_sum_same(const arma::mat& x) {
  const arma::uword dimension = x.n_rows;
  const arma::mat identity_mat = arma::eye(dimension, dimension);
  return arma::kron(x, identity_mat) + arma::kron(identity_mat, x);
}
//==============================================================================


//============================ kronecker(x, I) =================================
arma::mat kronecker_vector_identity(const arma::vec& x) {
  const arma::uword dimension = x.n_rows;
  return arma::kron(x, arma::eye(dimension, dimension));
}
//==============================================================================


//============================ kronecker(x, Y) =================================
arma::mat kronecker_vector_matrix(const arma::vec& x, const arma::mat& y) {
  return arma::kron(x, y);
}
//==============================================================================


//============================ X^{-1} y ========================================
arma::vec solve_chol_or_lu_vec(const arma::mat& X, const arma::vec& y) {
  arma::mat chol_upper;
  if (arma::chol(chol_upper, X)) {
    const arma::vec forward =
      arma::solve(arma::trimatl(chol_upper.t()), y);
    return arma::solve(arma::trimatu(chol_upper), forward);
  }
  arma::mat lu_lower, lu_upper, permutation_matrix;
  if (!arma::lu(lu_lower, lu_upper, permutation_matrix, X)) {
    Rcpp::stop("solve_chol_or_lu_vec: LU factorisation failed -- "
                 "matrix is singular or numerically unstable.");
  }
  const arma::vec permuted_rhs = permutation_matrix * y;
  arma::vec forward;
  if (!arma::solve(forward, arma::trimatl(lu_lower), permuted_rhs,
                   arma::solve_opts::no_approx)) {
    Rcpp::stop("solve_chol_or_lu_vec: LU forward solve failed -- "
                 "matrix is singular or numerically unstable.");
  }
  arma::vec ans;
  if (!arma::solve(ans, arma::trimatu(lu_upper), forward,
                   arma::solve_opts::no_approx)) {
    Rcpp::stop("solve_chol_or_lu_vec: LU back solve failed -- "
                 "matrix is singular or numerically unstable.");
  }
  return ans;
}
//==============================================================================


//============================ X^{-1} Y ========================================
arma::mat solve_chol_or_lu_mat(const arma::mat& X, const arma::mat& Y) {
  arma::mat chol_upper;
  if (arma::chol(chol_upper, X)) {
    const arma::mat forward =
      arma::solve(arma::trimatl(chol_upper.t()), Y);
    return arma::solve(arma::trimatu(chol_upper), forward);
  }
  arma::mat lu_lower, lu_upper, permutation_matrix;
  if (!arma::lu(lu_lower, lu_upper, permutation_matrix, X)) {
    Rcpp::stop("solve_chol_or_lu_mat: LU factorisation failed -- "
                 "matrix is singular or numerically unstable.");
  }
  const arma::mat permuted_rhs = permutation_matrix * Y;
  arma::mat forward;
  if (!arma::solve(forward, arma::trimatl(lu_lower), permuted_rhs,
                   arma::solve_opts::no_approx)) {
    Rcpp::stop("solve_chol_or_lu_mat: LU forward solve failed -- "
                 "matrix is singular or numerically unstable.");
  }
  arma::mat ans;
  if (!arma::solve(ans, arma::trimatu(lu_upper), forward,
                   arma::solve_opts::no_approx)) {
    Rcpp::stop("solve_chol_or_lu_mat: LU back solve failed -- "
                 "matrix is singular or numerically unstable.");
  }
  return ans;
}
//==============================================================================


//============================ symmetrize if close to symmetric ================
void symmetrize_if_close(arma::mat& A, double rel_tol) {
  const double scale = std::max(1.0, arma::abs(A).max());
  double max_asym = 0.0;
  const arma::uword n = A.n_rows;
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      const double d = std::abs(A(i, j) - A(j, i));
      if (d > max_asym) {
        max_asym = d;
      }
    }
  }
  if (max_asym > rel_tol * scale) {
    Rcpp::warning(
      "Fisher information matrix analogue has non-trivial asymmetry (%.2e); "
      "check for numerical issues.",
      max_asym / scale
    );
  }
  A = 0.5 * (A + A.t());
}
//==============================================================================


//============================ lambda from matrix blocks =======================
arma::vec lambda_from_blocks_chol_or_lu(const arma::mat& A,
                                        const arma::mat& lambda_matrix) {
  const arma::uword p = A.n_rows;
  arma::vec lambda_vector(p, arma::fill::zeros);
  arma::mat chol_upper;
  if (arma::chol(chol_upper, A)) {
    for (arma::uword r = 0; r < p; ++r) {
      const arma::mat block = lambda_matrix.rows(r * p, (r + 1) * p - 1);
      const arma::mat forward =
        arma::solve(arma::trimatl(chol_upper.t()), block);
      const arma::mat solved =
        arma::solve(arma::trimatu(chol_upper), forward);
      lambda_vector[r] = -0.5 * arma::trace(solved);
    }
    return lambda_vector;
  }
  arma::mat lu_lower, lu_upper, permutation_matrix;
  if (!arma::lu(lu_lower, lu_upper, permutation_matrix, A)) {
    Rcpp::stop("lambda_from_blocks_chol_or_lu: LU factorisation failed -- "
                 "matrix is singular or numerically unstable.");
  }
  for (arma::uword r = 0; r < p; ++r) {
    const arma::mat block = lambda_matrix.rows(r * p, (r + 1) * p - 1);
    const arma::mat permuted_block = permutation_matrix * block;
    arma::mat forward;
    if (!arma::solve(forward, arma::trimatl(lu_lower), permuted_block,
                     arma::solve_opts::no_approx)) {
      Rcpp::stop("lambda_from_blocks_chol_or_lu: LU forward solve failed -- "
                   "matrix is singular or numerically unstable.");
    }
    arma::mat solved;
    if (!arma::solve(solved, arma::trimatu(lu_upper), forward,
                     arma::solve_opts::no_approx)) {
      Rcpp::stop("lambda_from_blocks_chol_or_lu: LU back solve failed -- "
                   "matrix is singular or numerically unstable.");
    }
    lambda_vector[r] = -0.5 * arma::trace(solved);
  }
  return lambda_vector;
}
//==============================================================================
