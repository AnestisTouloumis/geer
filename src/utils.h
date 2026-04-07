#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

Rcpp::NumericVector arma2vec(const arma::vec& x);
arma::vec vec2arma(const Rcpp::NumericVector& x);
arma::mat subset_matrix(const arma::mat& x, const arma::vec& y);
arma::mat kappa_matrix(arma::uword dimension);
arma::mat kronecker_left_identity_kappa(const arma::mat& x);
arma::mat kronecker_identity_right_kappa(const arma::mat& x);
arma::mat kappa_right(const arma::mat& x);
arma::mat kronecker_sum_same(const arma::mat& x);
arma::mat kronecker_vector_identity(const arma::vec& x);
arma::mat kronecker_vector_matrix(const arma::vec& x, const arma::mat& y);
arma::vec solve_chol_or_lu_vec(const arma::mat& X, const arma::vec& y);
arma::mat solve_chol_or_lu_mat(const arma::mat& X, const arma::mat& Y);
arma::vec lambda_from_blocks_chol_or_lu(const arma::mat& A,
                                        const arma::mat& lambda_matrix);
inline void kron_self_matrix_into(arma::mat& kron_d_matrix_d_matrix_i,
                                  const arma::mat& d_matrix_i) {
  const arma::uword m = d_matrix_i.n_rows;
  const arma::uword p = d_matrix_i.n_cols;
  if (kron_d_matrix_d_matrix_i.n_rows != m * m ||
      kron_d_matrix_d_matrix_i.n_cols != p * p) {
    kron_d_matrix_d_matrix_i.set_size(m * m, p * p);
  }
  for (arma::uword j = 0; j < p; ++j) {
    const double* Dj = d_matrix_i.colptr(j);
    for (arma::uword k = 0; k < p; ++k) {
      const double* Dk = d_matrix_i.colptr(k);
      const arma::uword col_out = j * p + k;
      double* out = kron_d_matrix_d_matrix_i.colptr(col_out);

      for (arma::uword r = 0; r < m; ++r) {
        const double s = Dj[r];
        double* block = out + r * m;

        for (arma::uword i = 0; i < m; ++i) {
          block[i] = s * Dk[i];
        }
      }
    }
  }
}

#endif
