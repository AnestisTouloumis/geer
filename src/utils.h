#ifndef GEER_UTILS_H
#define GEER_UTILS_H

#include <RcppArmadillo.h>
#include <exception>
#include <functional>
#include <initializer_list>

struct GeerInputError : public Rcpp::exception {
  explicit GeerInputError(const char* message)
      : Rcpp::exception(message, false) {}
};

template <typename... Args>
[[noreturn]] inline void signal_input_error(const char* format, Args... args) {
  try {
    Rcpp::stop(format, args...);
  } catch (const Rcpp::exception& error) {
    throw GeerInputError(error.what());
  }
}

inline arma::uword upper_triangular_pairs(const arma::uword n) {
  return n * (n - 1) / 2;
}
inline arma::uword upper_triangular_pair_index(const arma::uword i,
                                               const arma::uword j,
                                               const arma::uword n) {
  return i * (2 * n - i - 1) / 2 + (j - i - 1);
}
arma::mat subset_matrix(const arma::mat& x, const arma::vec& y);
arma::mat vecdiag_matrix(const arma::uword dimension);
arma::mat kronecker_left_identity_vecdiag(const arma::mat& x);
arma::mat kronecker_identity_right_vecdiag(const arma::mat& x);
arma::mat vecdiag_right(const arma::mat& x);
inline arma::mat vecdiag_right_diag(const arma::vec& v) {
  const arma::uword dimension = v.n_elem;
  arma::mat result(dimension * dimension, dimension, arma::fill::zeros);
  for (arma::uword i = 0; i < dimension; ++i) {
    result(i * (dimension + 1), i) = v[i];
  }
  return result;
}
arma::mat kronecker_sum_same(const arma::mat& x);
class CholOrLuFactor {
 public:
  explicit CholOrLuFactor(const arma::mat& X,
                          const char* context = "solve_chol_or_lu",
                          bool allow_chol = true);
  arma::mat solve(const arma::mat& Y) const;
  arma::vec solve(const arma::vec& y) const;
private:
  template <typename T> T solve_impl(const T& rhs) const;
  const char* context_;
  bool use_chol_;
  bool symmetric_not_pd_;
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
[[noreturn]] void rethrow_with_cluster_context(const char* where,
                                               arma::uword cluster_index,
                                               double id_value,
                                               const std::exception& e);
arma::vec get_lambda_from_blocks(const arma::mat& A,
                                        const arma::mat& lambda_matrix);
void symmetrize_if_close(arma::mat& A, const char* context,
                         double rel_tol = 1e-10);
constexpr double leverage_tolerance = 1e-8;
void reset_lu_fallback_warning();
void check_working_covariance(const arma::mat& matrix, const char* context);
double get_max_cluster_leverage(const arma::mat& naive_root,
                                const arma::mat& v_inverse_d_matrix,
                                const arma::mat& d_matrix);
using vector_refs =
  std::initializer_list<std::reference_wrapper<const arma::vec>>;
void check_cluster_inputs(const arma::vec& id_vector,
                          vector_refs observation_vectors);
void check_cluster_inputs(const arma::vec& id_vector,
                          vector_refs observation_vectors,
                          const arma::mat& model_matrix);
inline void add_kronecker_self_t_s_d(arma::mat& out,
                                const arma::mat& d_matrix,
                                const arma::mat& s_matrix,
                                const double scale = 1.0) {
  const arma::uword cluster_size = d_matrix.n_rows;
  const arma::uword p = d_matrix.n_cols;
  const arma::mat sd_matrix = s_matrix * d_matrix;
  arma::mat p_matrix(p, p);
  for (arma::uword l = 0; l < p; ++l) {
    const arma::mat t_matrix =
      arma::reshape(sd_matrix.col(l), cluster_size, cluster_size);
    p_matrix = d_matrix.t() * t_matrix.t() * d_matrix;
    out.col(l) += scale * arma::vectorise(p_matrix.t());
  }
}

inline arma::vec kronecker_self_t_vec(const arma::mat& d_matrix,
                                 const arma::vec& w_vector) {
  const arma::uword cluster_size = d_matrix.n_rows;
  const arma::mat w_matrix =
    arma::reshape(w_vector, cluster_size, cluster_size);
  const arma::mat p_matrix = d_matrix.t() * w_matrix.t() * d_matrix;
  return arma::vectorise(p_matrix.t());
}

#endif
