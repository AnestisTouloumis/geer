#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>


namespace {
constexpr double symmetry_tolerance = 1e-10;
bool lu_fallback_warned = false;

double get_relative_asymmetry(const arma::mat& X) {
  if (!X.is_square()) {
    Rcpp::stop("get_relative_asymmetry failed: the matrix must be square");
  }
  if (X.is_empty()) return 0.0;
  if (!X.is_finite()) return NA_REAL;
  const double scale = arma::abs(X).max();
  if (scale == 0.0) return 0.0;
  const arma::mat scaled = X / scale;
  return arma::abs(scaled - scaled.t()).max();
}
}  // namespace


//============================ subset matrix x[y, y] ===========================
arma::mat subset_matrix(const arma::mat& x, const arma::vec& y) {
  const arma::uword subset_size = y.n_elem;
  const double upper = static_cast<double>(x.n_rows);
  std::vector<arma::uword> zero_based_indices(subset_size);
  for (arma::uword i = 0; i < subset_size; ++i) {
    if (!(y[i] >= 1.0 && y[i] <= upper)) {
      signal_input_error(
        "indices must lie in [1, %d]",
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


//============================ factorize once, solve many ======================
CholOrLuFactor::CholOrLuFactor(const arma::mat& X,
                               const char* context,
                               bool allow_chol)
    : context_(context), use_chol_(false), symmetric_not_pd_(false) {
  if (!X.is_finite()) {
    Rcpp::stop("%s failed: the matrix must be finite", context_);
  }
  if (allow_chol && X.is_square() && X.n_elem > 0 &&
      get_relative_asymmetry(X) <= symmetry_tolerance) {
    if (arma::chol(chol_upper_, X)) {
      use_chol_ = true;
      return;
    }
    symmetric_not_pd_ = true;
  }
  if (!arma::lu(lu_lower_, lu_upper_, permutation_, X)) {
    Rcpp::stop("%s failed: LU factorization failed; the matrix is "
      "singular or numerically unstable",
               context_);
  }
  if (!allow_chol) {
    const arma::vec abs_pivots = arma::abs(lu_upper_.diag());
    if (abs_pivots.min() <
        abs_pivots.max() * arma::datum::eps * static_cast<double>(X.n_rows)) {
      Rcpp::stop("%s failed: the matrix is numerically singular", context_);
    }
  }
}

template <typename T>
T CholOrLuFactor::solve_impl(const T& rhs) const {
  if (use_chol_) {
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
    Rcpp::stop("%s failed: LU forward solve failed; the matrix is "
      "singular or numerically unstable",
               context_);
  }
  T result;
  if (!arma::solve(result, arma::trimatu(lu_upper_), forward,
                   arma::solve_opts::no_approx)) {
    Rcpp::stop("%s failed: LU back solve failed; the matrix is "
      "singular or numerically unstable",
               context_);
  }
  if (symmetric_not_pd_ && !lu_fallback_warned) {
    lu_fallback_warned = true;
    Rcpp::Function warning_function("warning");
    warning_function(
      "a symmetric matrix that is not positive definite was solved by LU "
      "decomposition",
      Rcpp::Named("call.") = false);
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
  if (dynamic_cast<const GeerInputError*>(&e) != nullptr) {
    signal_input_error("%s; cluster %d (id = %g)", e.what(),
                       static_cast<int>(cluster_index + 1), id_value);
  }
  Rcpp::stop("%s failed: cluster %d (id = %g): %s", where,
             static_cast<int>(cluster_index + 1), id_value, e.what());
}
//==============================================================================


//============================ symmetrize if close to symmetric ================
void symmetrize_if_close(arma::mat& A, const char* context, double rel_tol) {
  const double relative_asymmetry = get_relative_asymmetry(A);
  if (ISNAN(relative_asymmetry)) {
    bool undefined = true;
    for (double value : A) {
      if (!ISNAN(value)) undefined = false;
    }
    if (undefined) return;
    Rcpp::stop("%s failed: the matrix must be finite", context);
  }
  if (relative_asymmetry > rel_tol) {
    Rcpp::stop("%s failed: the matrix has non-trivial asymmetry "
      "(%.2e)", context,
               relative_asymmetry);
  }
  if (!A.is_empty() &&
      arma::abs(A).max() > std::numeric_limits<double>::max() / 2.0) {
    A = 0.5 * A + 0.5 * A.t();
  } else {
    A = 0.5 * (A + A.t());
  }
}
//==============================================================================


//============================ lambda from matrix blocks =======================
arma::vec get_lambda_from_blocks(const arma::mat& A,
                                        const arma::mat& lambda_matrix) {
  const arma::uword p = A.n_rows;
  const CholOrLuFactor factor(A, "get_lambda_from_blocks");
  arma::vec lambda_vector(p, arma::fill::zeros);
  for (arma::uword r = 0; r < p; ++r) {
    const arma::mat block = lambda_matrix.rows(r * p, (r + 1) * p - 1);
    lambda_vector[r] = -0.5 * arma::trace(factor.solve(block));
  }
  return lambda_vector;
}
//==============================================================================


//============================ cluster inputs ==================================
void check_cluster_inputs(const arma::vec& id_vector,
                          vector_refs observation_vectors) {
  const arma::uword n = id_vector.n_elem;
  if (n == 0) {
    signal_input_error("input vectors must not be empty");
  }
  for (const arma::vec& observation_vector : observation_vectors) {
    if (observation_vector.n_elem != n) {
      signal_input_error("all observation-level inputs must have the same "
      "length");
    }
    if (!observation_vector.is_finite()) {
      signal_input_error("all observation-level inputs must be finite");
    }
  }
}


void check_cluster_inputs(const arma::vec& id_vector,
                          vector_refs observation_vectors,
                          const arma::mat& model_matrix) {
  check_cluster_inputs(id_vector, observation_vectors);
  if (model_matrix.n_rows != id_vector.n_elem) {
    signal_input_error("'model_matrix' must have one row per observation");
  }
  if (!model_matrix.is_finite()) {
    signal_input_error("'model_matrix' must be finite");
  }
}
//==============================================================================


//============================ maximum cluster leverage ========================
double get_max_cluster_leverage(const arma::mat& naive_root,
                                const arma::mat& v_inverse_d_matrix,
                                const arma::mat& d_matrix) {
  arma::mat cluster_leverage =
    naive_root * v_inverse_d_matrix.t() * d_matrix * naive_root.t();
  cluster_leverage = 0.5 * (cluster_leverage + cluster_leverage.t());
  return arma::max(arma::eig_sym(cluster_leverage));
}
//==============================================================================


//============================ reset LU fallback warning =======================
void reset_lu_fallback_warning() {
  lu_fallback_warned = false;
}
//==============================================================================


void check_working_covariance(const arma::mat& matrix, const char* context) {
  if (!matrix.is_finite()) {
    Rcpp::stop("%s failed: the working covariance matrix must be finite",
               context);
  }
  arma::mat root;
  if (!arma::chol(root, matrix)) {
    Rcpp::stop("%s failed: the working covariance matrix is not positive "
               "definite", context);
  }
}
