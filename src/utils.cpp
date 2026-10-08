#include "utils.h"
#include <algorithm>
#include <cmath>
#include <vector>


//============================ subset matrix x[y, y] ===========================
arma::mat subset_matrix(const arma::mat& x, const arma::vec& y) {
  const arma::uword subset_size = y.n_elem;
  const double upper = static_cast<double>(x.n_rows);
  std::vector<arma::uword> zero_based_indices(subset_size);
  for (arma::uword i = 0; i < subset_size; ++i) {
    // The negated comparison also rejects NaN.
    if (!(y[i] >= 1.0 && y[i] <= upper)) {
      Rcpp::stop(
        "subset_matrix: index out of range -- values must be in [1, %d].",
        static_cast<int>(x.n_rows)
      );
    }
    zero_based_indices[i] = static_cast<arma::uword>(y[i]) - 1;
  }
  arma::mat result(subset_size, subset_size);
  for (arma::uword j = 0; j < subset_size; ++j) {
    for (arma::uword i = 0; i < subset_size; ++i) {
      result(i, j) = x(zero_based_indices[i], zero_based_indices[j]);
    }
  }
  return result;
}
//==============================================================================


//============================ vecdiag matrix ==================================
arma::mat vecdiag_matrix(const arma::uword dimension) {
  arma::mat result(dimension * dimension, dimension, arma::fill::zeros);
  for (arma::uword j = 0; j < dimension; ++j) {
    result(j * dimension + j, j) = 1.0;
  }
  return result;
}
//==============================================================================


//============================ kronecker(X, I) * vecdiag =======================
arma::mat kronecker_left_identity_vecdiag(const arma::mat& x) {
  const arma::uword dimension = x.n_rows;
  arma::mat result(dimension * dimension, dimension, arma::fill::zeros);
  for (arma::uword j = 0; j < dimension; ++j) {
    const arma::uword base = j * dimension;
    for (arma::uword i = 0; i < dimension; ++i) {
      result(base + i, i) = x(j, i);
    }
  }
  return result;
}
//==============================================================================


//============================ kronecker(I, X) * vecdiag =======================
arma::mat kronecker_identity_right_vecdiag(const arma::mat& x) {
  const arma::uword dimension = x.n_rows;
  arma::mat result(dimension * dimension, dimension, arma::fill::zeros);
  for (arma::uword j = 0; j < dimension; ++j) {
    result.submat(j * dimension, j, (j + 1) * dimension - 1, j) = x.col(j);
  }
  return result;
}
//==============================================================================


//============================ vec-diag selector * X ===========================
arma::mat vecdiag_right(const arma::mat& x) {
  const arma::uword dimension = x.n_rows;
  arma::mat result(dimension * dimension, dimension, arma::fill::zeros);
  const arma::uvec rows =
    arma::regspace<arma::uvec>(0, dimension - 1) * (dimension + 1);
  result.rows(rows) = x;
  return result;
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


//============================ factorize once, solve many ======================
CholOrLuFactor::CholOrLuFactor(const arma::mat& X, const char* context)
  : context_(context), use_chol_(false) {
  if (arma::chol(chol_upper_, X)) {
    use_chol_ = true;
    return;
  }
  if (!arma::lu(lu_lower_, lu_upper_, permutation_, X)) {
    Rcpp::stop("%s: LU factorization failed -- "
                 "matrix is singular or numerically unstable.", context_);
  }
}

template <typename T>
T CholOrLuFactor::solve_impl(const T& rhs) const {
  if (use_chol_) {
    // Cholesky succeeded: X is positive-definite so triangular solves
    // cannot fail -- no explicit check needed.
    const T forward =
      arma::solve(arma::trimatl(chol_upper_.t()), rhs,
                  arma::solve_opts::no_approx);
    return arma::solve(arma::trimatu(chol_upper_), forward,
                       arma::solve_opts::no_approx);
  }
  const T permuted_rhs = permutation_ * rhs;
  T forward;
  if (!arma::solve(forward, arma::trimatl(lu_lower_), permuted_rhs,
                   arma::solve_opts::no_approx)) {
    Rcpp::stop("%s: LU forward solve failed -- "
                 "matrix is singular or numerically unstable.", context_);
  }
  T result;
  if (!arma::solve(result, arma::trimatu(lu_upper_), forward,
                   arma::solve_opts::no_approx)) {
    Rcpp::stop("%s: LU back solve failed -- "
                 "matrix is singular or numerically unstable.", context_);
  }
  return result;
}

arma::mat CholOrLuFactor::solve(const arma::mat& Y) const {
  return solve_impl<arma::mat>(Y);
}

arma::vec CholOrLuFactor::solve(const arma::vec& y) const {
  return solve_impl<arma::vec>(y);
}
//==============================================================================


//============================ X^{-1} y ========================================
arma::vec solve_chol_or_lu_vec(const arma::mat& X,
                               const arma::vec& y,
                               const char* context) {
  const CholOrLuFactor factor(X, context);
  return factor.solve(y);
}
//==============================================================================


//============================ X^{-1} Y ========================================
arma::mat solve_chol_or_lu_mat(const arma::mat& X,
                               const arma::mat& Y,
                               const char* context) {
  const CholOrLuFactor factor(X, context);
  return factor.solve(Y);
}
//==============================================================================


//============================ cluster error context ===========================
[[noreturn]] void rethrow_with_cluster_context(const char* where,
                                               arma::uword cluster_index,
                                               double id_value,
                                               const std::exception& e) {
  Rcpp::stop("%s: cluster %d (id = %g): %s",
             where,
             static_cast<int>(cluster_index + 1),
             id_value,
             e.what());
}
//==============================================================================


//============================ symmetrize if close to symmetric ================
void symmetrize_if_close(arma::mat& A, double rel_tol) {
  const double scale = std::max(1.0, arma::abs(A).max());
  const double max_asym = arma::abs(A - A.t()).max();
  if (max_asym > rel_tol * scale) {
    Rcpp::warning(
      "matrix has non-trivial asymmetry (%.2e) before symmetrization; "
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
  const CholOrLuFactor factor(A, "lambda_from_blocks_chol_or_lu");
  arma::vec lambda_vector(p, arma::fill::zeros);
  for (arma::uword r = 0; r < p; ++r) {
    const arma::mat block = lambda_matrix.rows(r * p, (r + 1) * p - 1);
    lambda_vector[r] = -0.5 * arma::trace(factor.solve(block));
  }
  return lambda_vector;
}
//==============================================================================
