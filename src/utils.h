#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>
#include <exception>

arma::mat subset_matrix(const arma::mat& x, const arma::vec& y);
arma::mat kappa_matrix(const arma::uword dimension);
arma::mat kronecker_left_identity_kappa(const arma::mat& x);
arma::mat kronecker_identity_right_kappa(const arma::mat& x);
arma::mat kappa_right(const arma::mat& x);
inline arma::mat kappa_right_diag(const arma::vec& v) {
  const arma::uword m = v.n_elem;
  arma::mat ans(m * m, m, arma::fill::zeros);
  for (arma::uword i = 0; i < m; ++i) {
    ans(i * (m + 1), i) = v[i];
  }
  return ans;
}
arma::mat kronecker_sum_same(const arma::mat& x);
arma::mat kronecker_vector_identity(const arma::vec& x);
// Factorizes a square matrix once (Cholesky when positive definite, LU with
// partial pivoting otherwise) so that several right-hand sides can be solved
// without repeating the factorization.
class CholOrLuFactor {
public:
  // `context` (a string literal) prefixes error messages so that a failed
  // factorization or solve names the calling routine.
  explicit CholOrLuFactor(const arma::mat& X,
                          const char* context = "solve_chol_or_lu");
  arma::mat solve(const arma::mat& Y) const;
  arma::vec solve(const arma::vec& y) const;
private:
  template <typename T> T solve_impl(const T& rhs) const;
  const char* context_;
  bool use_chol_;
  arma::mat chol_upper_;
  arma::mat lu_lower_;
  arma::mat lu_upper_;
  arma::mat permutation_;
};
arma::vec solve_chol_or_lu_vec(const arma::mat& X,
                               const arma::vec& y,
                               const char* context = "solve_chol_or_lu");
arma::mat solve_chol_or_lu_mat(const arma::mat& X,
                               const arma::mat& Y,
                               const char* context = "solve_chol_or_lu");
// Re-throws `e` as an R error that names the routine and the cluster
// (1-based position and id value) in which it occurred.
[[noreturn]] void rethrow_with_cluster_context(const char* where,
                                               arma::uword cluster_index,
                                               double id_value,
                                               const std::exception& e);
arma::vec lambda_from_blocks_chol_or_lu(const arma::mat& A,
                                        const arma::mat& lambda_matrix);
void symmetrize_if_close(arma::mat& A, double rel_tol = 1e-10);
// out += scale * kron(D, D)' * S * D, where S is (m^2 x m) and out is
// (p^2 x p), without forming kron(D, D). Row j * p + k of out matches the
// (j, k) pair of the Kronecker product.
inline void add_kron_self_t_s_d(arma::mat& out,
                                const arma::mat& d_matrix,
                                const arma::mat& s_matrix,
                                const double scale = 1.0) {
  const arma::uword m = d_matrix.n_rows;
  const arma::uword p = d_matrix.n_cols;
  const arma::mat sd_matrix = s_matrix * d_matrix;
  arma::mat p_matrix(p, p);
  for (arma::uword l = 0; l < p; ++l) {
    const arma::mat t_matrix(const_cast<double*>(sd_matrix.colptr(l)), m, m,
                             false, true);
    p_matrix = d_matrix.t() * t_matrix.t() * d_matrix;
    double* out_col = out.colptr(l);
    for (arma::uword j = 0; j < p; ++j) {
      for (arma::uword k = 0; k < p; ++k) {
        out_col[j * p + k] += scale * p_matrix(j, k);
      }
    }
  }
}

// Returns kron(D, D)' * w for a vector w of length m^2, without forming
// kron(D, D). Entry j * p + k matches the (j, k) pair of the Kronecker product.
inline arma::vec kron_self_t_vec(const arma::mat& d_matrix,
                                 const arma::vec& w_vector) {
  const arma::uword m = d_matrix.n_rows;
  const arma::uword p = d_matrix.n_cols;
  const arma::mat w_matrix(const_cast<double*>(w_vector.memptr()), m, m,
                           false, true);
  const arma::mat p_matrix = d_matrix.t() * w_matrix.t() * d_matrix;
  arma::vec ans(p * p);
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = 0; k < p; ++k) {
      ans[j * p + k] = p_matrix(j, k);
    }
  }
  return ans;
}

#endif
