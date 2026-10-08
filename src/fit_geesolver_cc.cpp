#include "nuisance_estimation_cc.h"
#include "working_covariance_cc.h"
#include "cluster_utils.h"
#include "utils.h"
#include "link_functions.h"
#include "variance_functions.h"
#include "covariance_matrices.h"
#include "method_codes.h"
#include <algorithm>
#include <cmath>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>


namespace {

//============================ update beta - gee ===============================
arma::vec update_beta_gee_cc(const arma::vec& y_vector,
                             const arma::mat& model_matrix,
                             const arma::vec& id_vector,
                             const std::vector<Cluster>& clusters,
                             const arma::vec& repeated_vector,
                             const arma::vec& weights_vector,
                             const char* link,
                             const char* family,
                             const arma::vec& beta_vector,
                             const arma::vec& mu_vector,
                             const arma::vec& eta_vector,
                             const char* correlation_structure,
                             const arma::vec& alpha_vector,
                             const double phi) {
  const LinkCode link_code = parse_link(link);
  const FamilyCode family_code = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
  const arma::vec delta_vector = geer::link_derivative_1(link_code, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  arma::mat d_matrix_i;
  arma::mat v_matrix_i;
  arma::mat v_matrix_inverse_d_matrix_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      const auto delta_vector_i = delta_vector.subvec(first_row, last_row);
      const auto s_vector_i = s_vector.subvec(first_row, last_row);
      d_matrix_i = model_matrix.rows(first_row, last_row);
      d_matrix_i.each_col() %= delta_vector_i;
      v_matrix_i = get_v_matrix_cc(family_code,
                                   mu_vector.subvec(first_row, last_row),
                                   repeated_vector.subvec(first_row, last_row),
                                   phi,
                                   correlation_matrix,
                                   weights_vector.subvec(first_row, last_row));
      const CholOrLuFactor v_factor(v_matrix_i);
      v_matrix_inverse_d_matrix_i =
        v_factor.solve(d_matrix_i);
      d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
      information_matrix += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector_i;
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("update_beta_gee_cc", cluster_index, id_vector[cluster.start], e);
    }
  }
  symmetrize_if_close(information_matrix, 1e-10);
  const arma::vec update_step =
    solve_chol_or_lu_vec(information_matrix, u_vector,
        "update_beta_gee_cc: naive information matrix");
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - naive =============================
arma::vec update_beta_naive_cc(const arma::vec& y_vector,
                               const arma::mat& model_matrix,
                               const arma::vec& id_vector,
                               const std::vector<Cluster>& clusters,
                               const arma::vec& repeated_vector,
                               const arma::vec& weights_vector,
                               const char* link,
                               const char* family,
                               const arma::vec& beta_vector,
                               const arma::vec& mu_vector,
                               const arma::vec& eta_vector,
                               const char* correlation_structure,
                               const arma::vec& alpha_vector,
                               const double phi) {
  const LinkCode link_code = parse_link(link);
  const FamilyCode family_code = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat lambda_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = geer::link_derivative_1(link_code, eta_vector);
  const arma::vec delta_star_vector =
    geer::link_derivative_2(link_code, eta_vector) / arma::square(delta_vector);
  const arma::vec variance_vector = geer::variance_function(family_code, mu_vector);
  const arma::vec alpha_star_vector =
    -0.5 * geer::variance_derivative_1(family_code, mu_vector) / variance_vector;
  const arma::vec alpha_star_plus_delta_star_vector =
    alpha_star_vector + delta_star_vector;
  const arma::vec s_vector = y_vector - mu_vector;
  const arma::mat correlation_matrix =
    get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
  arma::mat d_matrix_i;
  arma::mat v_matrix_i;
  arma::mat v_matrix_inverse_d_matrix_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::mat v_matrix_inverse_alpha_plus_delta_star_diag_i;
  arma::mat v_matrix_inverse_delta_star_diag_i;
  arma::mat alpha_star_plus_delta_star_matrix_i;
  arma::mat delta_star_matrix_i;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.start;
      const arma::uword last_row = cluster.end - 1;
      const arma::uword cluster_size = cluster.end - cluster.start;
      const auto delta_vector_i = delta_vector.subvec(first_row, last_row);
      const auto s_vector_i = s_vector.subvec(first_row, last_row);
      const auto delta_star_vector_i =
        delta_star_vector.subvec(first_row, last_row);
      const auto alpha_star_plus_delta_star_vector_i =
        alpha_star_plus_delta_star_vector.subvec(first_row, last_row);
      d_matrix_i = model_matrix.rows(first_row, last_row);
      d_matrix_i.each_col() %= delta_vector_i;
      v_matrix_i = get_v_matrix_cc(family_code,
                                   mu_vector.subvec(first_row, last_row),
                                   repeated_vector.subvec(first_row, last_row),
                                   phi,
                                   correlation_matrix,
                                   weights_vector.subvec(first_row, last_row));
      const CholOrLuFactor v_factor(v_matrix_i);

      v_matrix_inverse_d_matrix_i =
        v_factor.solve(d_matrix_i);
      d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
      information_matrix += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector_i;
      if (alpha_star_plus_delta_star_matrix_i.n_rows != cluster_size ||
          alpha_star_plus_delta_star_matrix_i.n_cols != cluster_size) {
        alpha_star_plus_delta_star_matrix_i.set_size(cluster_size, cluster_size);
      }
      alpha_star_plus_delta_star_matrix_i.zeros();
      alpha_star_plus_delta_star_matrix_i.diag() =
        alpha_star_plus_delta_star_vector_i;
      if (delta_star_matrix_i.n_rows != cluster_size || delta_star_matrix_i.n_cols != cluster_size) {
        delta_star_matrix_i.set_size(cluster_size, cluster_size);
      }
      delta_star_matrix_i.zeros();
      delta_star_matrix_i.diag() = delta_star_vector_i;
      v_matrix_inverse_alpha_plus_delta_star_diag_i =
        v_factor.solve(alpha_star_plus_delta_star_matrix_i);
      v_matrix_inverse_delta_star_diag_i =
        v_factor.solve(delta_star_matrix_i);
      add_kron_self_t_s_d(
        lambda_matrix,
        d_matrix_i,
        vecdiag_right(v_matrix_inverse_alpha_plus_delta_star_diag_i.t()) -
          kronecker_identity_right_vecdiag(v_matrix_inverse_alpha_plus_delta_star_diag_i) -
          kronecker_left_identity_vecdiag(v_matrix_inverse_delta_star_diag_i));
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("update_beta_naive_cc", cluster_index, id_vector[cluster.start], e);
    }
  }
  symmetrize_if_close(information_matrix, 1e-10);
  const arma::vec lambda_vector =
    lambda_from_blocks_chol_or_lu(information_matrix, lambda_matrix);
  const arma::vec update_step =
    solve_chol_or_lu_vec(information_matrix, u_vector + lambda_vector,
      "update_beta_naive_cc: naive information matrix");
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - robust  ===========================
arma::vec update_beta_robust_cc(const arma::vec& y_vector,
                                const arma::mat& model_matrix,
                                const arma::vec& id_vector,
                                const std::vector<Cluster>& clusters,
                                const arma::vec& repeated_vector,
                                const arma::vec& weights_vector,
                                const char* link,
                                const char* family,
                                const arma::vec& beta_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& eta_vector,
                                const char* correlation_structure,
                                const arma::vec& alpha_vector,
                                const double phi) {
  const LinkCode link_code = parse_link(link);
  const FamilyCode family_code = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat partial_derivatives_matrix(params_no * params_no,
                                       params_no,
                                       arma::fill::zeros);
  arma::mat second_derivatives_matrix(params_no * params_no,
                                      params_no,
                                      arma::fill::zeros);
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = geer::link_derivative_1(link_code, eta_vector);
  const arma::vec delta_star_vector =
    geer::link_derivative_2(link_code, eta_vector) / arma::square(delta_vector);
  const arma::vec variance_vector = geer::variance_function(family_code, mu_vector);
  const arma::vec alpha_star_vector =
    -0.5 * geer::variance_derivative_1(family_code, mu_vector) / variance_vector;
    const arma::vec alpha_star_plus_delta_star_vector =
    alpha_star_vector + delta_star_vector;
    const arma::vec s_vector = y_vector - mu_vector;
    const arma::mat correlation_matrix =
      get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
    arma::mat d_matrix_i;
    arma::mat v_matrix_i;
    arma::mat v_matrix_inverse_d_matrix_i;
    arma::mat d_matrix_trans_v_matrix_inverse_i;
    arma::vec u_vector_i(params_no, arma::fill::zeros);
    arma::mat alpha_star_plus_delta_star_matrix_i;
    arma::mat alpha_star_matrix_i;
    arma::mat v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i;
    arma::mat v_matrix_inverse_alpha_star_matrix_i;
    for (arma::uword cluster_index = 0;
         cluster_index < clusters.size();
         ++cluster_index) {
      const auto& cluster = clusters[cluster_index];
      try {
        const arma::uword first_row = cluster.start;
        const arma::uword last_row = cluster.end - 1;
        const arma::uword cluster_size = cluster.end - cluster.start;
        const auto delta_vector_i = delta_vector.subvec(first_row, last_row);
        const auto s_vector_i = s_vector.subvec(first_row, last_row);
        const auto alpha_star_vector_i =
          alpha_star_vector.subvec(first_row, last_row);
        const auto alpha_star_plus_delta_star_vector_i =
          alpha_star_plus_delta_star_vector.subvec(first_row, last_row);
        d_matrix_i = model_matrix.rows(first_row, last_row);
        d_matrix_i.each_col() %= delta_vector_i;
        v_matrix_i = get_v_matrix_cc(family_code,
                                     mu_vector.subvec(first_row, last_row),
                                     repeated_vector.subvec(first_row, last_row),
                                     phi,
                                     correlation_matrix,
                                     weights_vector.subvec(first_row, last_row));
        const CholOrLuFactor v_factor(v_matrix_i);
        v_matrix_inverse_d_matrix_i =
          v_factor.solve(d_matrix_i);
        d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
        information_matrix += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
        u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector_i;
        u_vector += u_vector_i;
        meat_matrix += u_vector_i * u_vector_i.t();
        if (alpha_star_plus_delta_star_matrix_i.n_rows != cluster_size ||
            alpha_star_plus_delta_star_matrix_i.n_cols != cluster_size) {
          alpha_star_plus_delta_star_matrix_i.set_size(cluster_size, cluster_size);
        }
        alpha_star_plus_delta_star_matrix_i.zeros();
        alpha_star_plus_delta_star_matrix_i.diag() =
          alpha_star_plus_delta_star_vector_i;
        if (alpha_star_matrix_i.n_rows != cluster_size || alpha_star_matrix_i.n_cols != cluster_size) {
          alpha_star_matrix_i.set_size(cluster_size, cluster_size);
        }
        alpha_star_matrix_i.zeros();
        alpha_star_matrix_i.diag() = alpha_star_vector_i;
        v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i =
          v_factor.solve(alpha_star_plus_delta_star_matrix_i);
        v_matrix_inverse_alpha_star_matrix_i =
          v_factor.solve(alpha_star_matrix_i);
        const arma::mat
        vecdiag_matrix_delta_star_matrix_plus_alpha_star_matrix_v_matrix_inverse_i =
          vecdiag_right(v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i.t());
        add_kron_self_t_s_d(
          second_derivatives_matrix,
          d_matrix_i,
          vecdiag_matrix_delta_star_matrix_plus_alpha_star_matrix_v_matrix_inverse_i +
            kronecker_left_identity_vecdiag(
              v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i +
                v_matrix_inverse_alpha_star_matrix_i
            ) +
            kronecker_identity_right_vecdiag(
              v_matrix_inverse_alpha_star_matrix_plus_delta_star_matrix_i
            ),
          -1.0);
        partial_derivatives_matrix +=
          kron_self_t_vec(
            d_matrix_i,
            (vecdiag_matrix_delta_star_matrix_plus_alpha_star_matrix_v_matrix_inverse_i +
              kronecker_left_identity_vecdiag(v_matrix_inverse_alpha_star_matrix_i)) *
              s_vector_i) *
          u_vector_i.t();
      } catch (const std::exception& e) {
        rethrow_with_cluster_context("update_beta_robust_cc", cluster_index, id_vector[cluster.start], e);
      }
    }
    symmetrize_if_close(information_matrix, 1e-10);
    const arma::mat naive_covariance_meat_matrix =
      solve_chol_or_lu_mat(information_matrix, meat_matrix,
        "update_beta_robust_cc: naive information matrix");
    const arma::mat robust_matrix =
      solve_chol_or_lu_mat(information_matrix, naive_covariance_meat_matrix.t(),
        "update_beta_robust_cc: naive information matrix");
    arma::vec lambda_vector(params_no, arma::fill::zeros);
    for (arma::uword r = 0; r < params_no; ++r) {
      const arma::mat first_block =
        partial_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
      const arma::mat second_block =
        second_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
      lambda_vector[r] =
        -(arma::trace(solve_chol_or_lu_mat(information_matrix, first_block,
        "update_beta_robust_cc: naive information matrix")) +
        0.5 * arma::trace(robust_matrix * second_block));
    }
    return beta_vector +
      solve_chol_or_lu_vec(information_matrix, u_vector + lambda_vector,
        "update_beta_robust_cc: naive information matrix");
}
//==============================================================================


//============================ update beta - empirical =========================
arma::vec update_beta_empirical_cc(const arma::vec& y_vector,
                                   const arma::mat& model_matrix,
                                   const arma::vec& id_vector,
                                   const std::vector<Cluster>& clusters,
                                   const arma::vec& repeated_vector,
                                   const arma::vec& weights_vector,
                                   const char* link,
                                   const char* family,
                                   const arma::vec& beta_vector,
                                   const arma::vec& mu_vector,
                                   const arma::vec& eta_vector,
                                   const char* correlation_structure,
                                   const arma::vec& alpha_vector,
                                   const double phi) {
  const LinkCode link_code = parse_link(link);
  const FamilyCode family_code = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat partial_derivatives_matrix(params_no * params_no,
                                       params_no,
                                       arma::fill::zeros);
  arma::mat second_derivatives_matrix(params_no * params_no,
                                      params_no,
                                      arma::fill::zeros);
  arma::mat observed_fisher_info_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec s_vector = y_vector - mu_vector;
  const arma::vec delta_vector = geer::link_derivative_1(link_code, eta_vector);
  const arma::vec mueta2_vector = geer::link_derivative_2(link_code, eta_vector);
  const arma::vec delta_star_vector =
    mueta2_vector / arma::square(delta_vector);
  const arma::vec delta_tilde_star_vector =
    (delta_vector % geer::link_derivative_3(link_code, eta_vector) -
    2.0 * arma::square(mueta2_vector)) /
      arma::square(arma::square(delta_vector));
  const arma::vec variance_vector = geer::variance_function(family_code, mu_vector);
  const arma::vec variancemu_vector = geer::variance_derivative_1(family_code, mu_vector);
  const arma::vec alpha_star_vector =
    -0.5 * variancemu_vector / variance_vector;
    const arma::vec alpha_tilde_star_vector =
    0.5 * (arma::square(variancemu_vector) / variance_vector -
    geer::variance_derivative_2(family_code, mu_vector)) / variance_vector;
    const arma::vec alpha_star_plus_delta_star_vector =
      alpha_star_vector + delta_star_vector;
    const arma::mat correlation_matrix =
      get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
    arma::mat d_matrix_i;
    arma::mat v_matrix_i;
    arma::mat v_matrix_inverse_d_matrix_i;
    arma::mat v_matrix_inverse_i;
    arma::mat epsilon_matrix_i;
    arma::mat observed_fisher_info_matrix_i;
    arma::mat alpha_star_plus_delta_star_matrix_i;
    arma::mat temp_diagonal_matrix;
    arma::mat hessian_term2;
    arma::mat hessian_term3;
    arma::mat diagonal_correction;
    arma::vec u_vector_i(params_no, arma::fill::zeros);
    arma::vec weighted_residuals_i;
    for (arma::uword cluster_index = 0;
         cluster_index < clusters.size();
         ++cluster_index) {
      const auto& cluster = clusters[cluster_index];
      try {
        const arma::uword first_row = cluster.start;
        const arma::uword last_row = cluster.end - 1;
        const arma::uword cluster_size = cluster.end - cluster.start;
        const arma::vec s_vector_i = s_vector.subvec(first_row, last_row);
        const arma::vec alpha_star_vector_i = alpha_star_vector.subvec(first_row, last_row);
        const arma::vec alpha_tilde_star_vector_i =
          alpha_tilde_star_vector.subvec(first_row, last_row);
        const arma::vec alpha_star_plus_delta_star_vector_i =
          alpha_star_plus_delta_star_vector.subvec(first_row, last_row);
        d_matrix_i = model_matrix.rows(first_row, last_row);
        d_matrix_i.each_col() %= delta_vector.subvec(first_row, last_row);

        v_matrix_i = get_v_matrix_cc(family_code,
                                     mu_vector.subvec(first_row, last_row),
                                     repeated_vector.subvec(first_row, last_row),
                                     phi,
                                     correlation_matrix,
                                     weights_vector.subvec(first_row, last_row));
        const CholOrLuFactor v_factor(v_matrix_i);
        v_matrix_inverse_d_matrix_i =
          v_factor.solve(d_matrix_i);
        u_vector_i = v_matrix_inverse_d_matrix_i.t() * s_vector_i;
        u_vector += u_vector_i;
        meat_matrix += u_vector_i * u_vector_i.t();
        v_matrix_inverse_i =
          v_factor.solve(arma::mat(arma::eye(cluster_size, cluster_size)));
        weighted_residuals_i = v_matrix_inverse_i * s_vector_i;
        epsilon_matrix_i.set_size(cluster_size, cluster_size);
        epsilon_matrix_i.zeros();
        epsilon_matrix_i.diag() =
          weighted_residuals_i % alpha_star_plus_delta_star_vector_i;
        temp_diagonal_matrix.set_size(cluster_size, cluster_size);
        temp_diagonal_matrix.zeros();
        temp_diagonal_matrix.diag() = s_vector_i % alpha_star_vector_i - 1.0;
        epsilon_matrix_i += v_matrix_inverse_i * temp_diagonal_matrix;
        observed_fisher_info_matrix_i = d_matrix_i.t() * epsilon_matrix_i * d_matrix_i;
        observed_fisher_info_matrix -= observed_fisher_info_matrix_i;
        partial_derivatives_matrix +=
          arma::vectorise(observed_fisher_info_matrix_i.t()) * u_vector_i.t();
        alpha_star_plus_delta_star_matrix_i.set_size(cluster_size, cluster_size);
        alpha_star_plus_delta_star_matrix_i.zeros();
        alpha_star_plus_delta_star_matrix_i.diag() =
          alpha_star_plus_delta_star_vector_i;
        const arma::mat hessian_term1 =
          alpha_star_plus_delta_star_matrix_i * epsilon_matrix_i;
        const arma::vec hessian_diag2 =
          alpha_star_plus_delta_star_vector_i %
          (alpha_star_plus_delta_star_vector_i + alpha_star_vector_i) %
          weighted_residuals_i;
        hessian_term2.set_size(cluster_size, cluster_size);
        hessian_term2.zeros();
        hessian_term2.diag() = -hessian_diag2;
        const arma::vec hessian_diag3 =
          (alpha_tilde_star_vector_i + delta_tilde_star_vector.subvec(first_row, last_row)) %
          weighted_residuals_i;
        hessian_term3.set_size(cluster_size, cluster_size);
        hessian_term3.zeros();
        hessian_term3.diag() = hessian_diag3;
        const arma::mat vecdiag_correction =
          vecdiag_right(hessian_term1 + hessian_term2 + hessian_term3);
        const arma::mat right_kronecker_correction =
          kronecker_identity_right_vecdiag(
            epsilon_matrix_i.t() * alpha_star_plus_delta_star_matrix_i
          );
        diagonal_correction.set_size(cluster_size, cluster_size);
        diagonal_correction.zeros();
        diagonal_correction.diag() =
          s_vector_i % alpha_tilde_star_vector_i - alpha_star_vector_i;
        const arma::mat left_kronecker_correction =
          kronecker_left_identity_vecdiag(
            epsilon_matrix_i * alpha_star_plus_delta_star_matrix_i +
              v_matrix_inverse_i * diagonal_correction
          );
        add_kron_self_t_s_d(second_derivatives_matrix,
                            d_matrix_i,
                            vecdiag_correction +
                              right_kronecker_correction +
                              left_kronecker_correction);
      } catch (const std::exception& e) {
        rethrow_with_cluster_context("update_beta_empirical_cc", cluster_index, id_vector[cluster.start], e);
      }
    }
    arma::mat robust_matrix(params_no, params_no, arma::fill::zeros);
    arma::vec lambda_vector(params_no, arma::fill::zeros);
    arma::mat lu_lower, lu_upper, permutation_matrix;
    // arma::lu() almost always returns true even for singular matrices; the
    // real singularity signal is a near-zero pivot on the diagonal of U.
    if (!arma::lu(lu_lower, lu_upper, permutation_matrix,
                  observed_fisher_info_matrix)) {
      Rcpp::stop("update_beta_empirical_cc: LU factorization failed -- "
                   "observed Fisher information is singular or numerically unstable.");
    }
    const arma::vec abs_pivots = arma::abs(lu_upper.diag());
    if (abs_pivots.min() <
      abs_pivots.max() * arma::datum::eps * static_cast<double>(params_no)) {
      Rcpp::stop("update_beta_empirical_cc: observed Fisher information matrix "
                   "is numerically singular -- check for collinearity or near-zero "
                   "cluster variances.");
    }
    auto solve_observed_fisher_matrix = [&](const arma::mat& rhs_matrix) -> arma::mat {
      const arma::mat permuted_rhs = permutation_matrix * rhs_matrix;
      arma::mat lu_forward;
      if (!arma::solve(lu_forward,
                       arma::trimatl(lu_lower),
                       permuted_rhs,
                       arma::solve_opts::no_approx)) {
        Rcpp::stop("update_beta_empirical_cc: forward LU solve failed.");
      }
      arma::mat result;
      if (!arma::solve(result,
                       arma::trimatu(lu_upper),
                       lu_forward,
                       arma::solve_opts::no_approx)) {
        Rcpp::stop("update_beta_empirical_cc: backward LU solve failed.");
      }
      return result;
    };
    auto solve_observed_fisher_vector = [&](const arma::vec& rhs_vector) -> arma::vec {
      const arma::vec permuted_rhs = permutation_matrix * rhs_vector;
      arma::vec lu_forward;
      if (!arma::solve(lu_forward,
                       arma::trimatl(lu_lower),
                       permuted_rhs,
                       arma::solve_opts::no_approx)) {
        Rcpp::stop("update_beta_empirical_cc: forward LU solve failed.");
      }
      arma::vec result;
      if (!arma::solve(result,
                       arma::trimatu(lu_upper),
                       lu_forward,
                       arma::solve_opts::no_approx)) {
        Rcpp::stop("update_beta_empirical_cc: backward LU solve failed.");
      }
      return result;
    };
    const arma::mat fisher_inv_meat = solve_observed_fisher_matrix(meat_matrix);
    robust_matrix = solve_observed_fisher_matrix(fisher_inv_meat.t());
    for (arma::uword param_idx = 0; param_idx < params_no; ++param_idx) {
      const arma::uword block_start = param_idx * params_no;
      const arma::uword block_end = (param_idx + 1) * params_no - 1;
      const arma::mat partial_deriv_block =
        partial_derivatives_matrix.rows(block_start, block_end);
      const arma::mat second_deriv_block =
        second_derivatives_matrix.rows(block_start, block_end);
      lambda_vector[param_idx] =
        -(arma::trace(solve_observed_fisher_matrix(partial_deriv_block)) +
        0.5 * arma::trace(robust_matrix * second_deriv_block));
    }
    const arma::vec update_step =
      solve_observed_fisher_vector(u_vector + lambda_vector);
    return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - jeffreys ==========================
arma::vec update_beta_jeffreys_cc(const arma::vec& y_vector,
                                  const arma::mat& model_matrix,
                                  const arma::vec& id_vector,
                                  const std::vector<Cluster>& clusters,
                                  const arma::vec& repeated_vector,
                                  const arma::vec& weights_vector,
                                  const char* link,
                                  const char* family,
                                  const arma::vec& beta_vector,
                                  const arma::vec& mu_vector,
                                  const arma::vec& eta_vector,
                                  const char* correlation_structure,
                                  const arma::vec& alpha_vector,
                                  const double phi,
                                  const double jeffreys_power) {
  const LinkCode link_code = parse_link(link);
  const FamilyCode family_code = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max =
    static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat lambda_matrix(params_no * params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = geer::link_derivative_1(link_code, eta_vector);
  const arma::vec delta_star_vector =
    geer::link_derivative_2(link_code, eta_vector) / arma::square(delta_vector);
  const arma::vec variance_vector = geer::variance_function(family_code, mu_vector);
  const arma::vec alpha_star_vector =
    -0.5 * geer::variance_derivative_1(family_code, mu_vector) / variance_vector;
    const arma::vec alpha_star_plus_delta_star_vector =
    alpha_star_vector + delta_star_vector;
    const arma::vec s_vector = y_vector - mu_vector;
    const arma::mat correlation_matrix =
      get_correlation_matrix(correlation_structure, alpha_vector, repeated_max);
    arma::mat d_matrix_i;
    arma::mat v_matrix_i;
    arma::mat v_matrix_inverse_d_matrix_i;
    arma::mat d_matrix_trans_v_matrix_inverse_i;
    arma::mat alpha_star_plus_delta_star_matrix_i;
    arma::mat v_matrix_inverse_alpha_star_plus_delta_star_matrix_i;
    for (arma::uword cluster_index = 0;
         cluster_index < clusters.size();
         ++cluster_index) {
      const auto& cluster = clusters[cluster_index];
      try {
        const arma::uword first_row = cluster.start;
        const arma::uword last_row = cluster.end - 1;
        const arma::uword cluster_size = cluster.end - cluster.start;
        const auto delta_vector_i = delta_vector.subvec(first_row, last_row);
        const auto s_vector_i = s_vector.subvec(first_row, last_row);
        const auto alpha_star_plus_delta_star_vector_i =
          alpha_star_plus_delta_star_vector.subvec(first_row, last_row);

        d_matrix_i = model_matrix.rows(first_row, last_row);
        d_matrix_i.each_col() %= delta_vector_i;
        v_matrix_i = get_v_matrix_cc(family_code,
                                     mu_vector.subvec(first_row, last_row),
                                     repeated_vector.subvec(first_row, last_row),
                                     phi,
                                     correlation_matrix,
                                     weights_vector.subvec(first_row, last_row));
        const CholOrLuFactor v_factor(v_matrix_i);
        v_matrix_inverse_d_matrix_i =
          v_factor.solve(d_matrix_i);
        d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
        information_matrix += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
        u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector_i;
        if (alpha_star_plus_delta_star_matrix_i.n_rows != cluster_size ||
            alpha_star_plus_delta_star_matrix_i.n_cols != cluster_size) {
          alpha_star_plus_delta_star_matrix_i.set_size(cluster_size, cluster_size);
        }
        alpha_star_plus_delta_star_matrix_i.zeros();
        alpha_star_plus_delta_star_matrix_i.diag() =
          alpha_star_plus_delta_star_vector_i;
        v_matrix_inverse_alpha_star_plus_delta_star_matrix_i =
          v_factor.solve(alpha_star_plus_delta_star_matrix_i);
        add_kron_self_t_s_d(
          lambda_matrix,
          d_matrix_i,
          kronecker_left_identity_vecdiag(
            v_matrix_inverse_alpha_star_plus_delta_star_matrix_i
          ) +
            kronecker_identity_right_vecdiag(
              v_matrix_inverse_alpha_star_plus_delta_star_matrix_i
            ));
      } catch (const std::exception& e) {
        rethrow_with_cluster_context("update_beta_jeffreys_cc", cluster_index, id_vector[cluster.start], e);
      }
    }
    symmetrize_if_close(information_matrix, 1e-10);
    const arma::mat naive_covariance =
      solve_chol_or_lu_mat(information_matrix,
                           arma::eye(params_no, params_no),
        "update_beta_jeffreys_cc: naive information matrix");
    const arma::vec naive_covariance_vectorized = arma::vectorise(naive_covariance);
    const arma::vec lambda_vector =
      jeffreys_power * (lambda_matrix.t() * naive_covariance_vectorized);
    return beta_vector +
      solve_chol_or_lu_vec(information_matrix, u_vector + lambda_vector,
        "update_beta_jeffreys_cc: naive information matrix");
}
//==============================================================================


//============================ update beta =====================================
arma::vec update_beta_cc(const arma::vec& y_vector,
                         const arma::mat& model_matrix,
                         const arma::vec& id_vector,
                         const std::vector<Cluster>& clusters,
                         const arma::vec& repeated_vector,
                         const arma::vec& weights_vector,
                         const char* link,
                         const char* family,
                         const arma::vec& beta_vector,
                         const arma::vec& mu_vector,
                         const arma::vec& eta_vector,
                         const char* correlation_structure,
                         const arma::vec& alpha_vector,
                         const double phi,
                         const double jeffreys_power,
                         const char* method) {
  switch (method_code(method)) {
  case MethodCode::gee:
    return update_beta_gee_cc(y_vector,
                              model_matrix,
                              id_vector,
                              clusters,
                              repeated_vector,
                              weights_vector,
                              link,
                              family,
                              beta_vector,
                              mu_vector,
                              eta_vector,
                              correlation_structure,
                              alpha_vector,
                              phi);
  case MethodCode::brgee_naive:
    return update_beta_naive_cc(y_vector,
                                model_matrix,
                                id_vector,
                                clusters,
                                repeated_vector,
                                weights_vector,
                                link,
                                family,
                                beta_vector,
                                mu_vector,
                                eta_vector,
                                correlation_structure,
                                alpha_vector,
                                phi);
  case MethodCode::brgee_robust:
    return update_beta_robust_cc(y_vector,
                                 model_matrix,
                                 id_vector,
                                 clusters,
                                 repeated_vector,
                                 weights_vector,
                                 link,
                                 family,
                                 beta_vector,
                                 mu_vector,
                                 eta_vector,
                                 correlation_structure,
                                 alpha_vector,
                                 phi);
  case MethodCode::brgee_empirical:
    return update_beta_empirical_cc(y_vector,
                                    model_matrix,
                                    id_vector,
                                    clusters,
                                    repeated_vector,
                                    weights_vector,
                                    link,
                                    family,
                                    beta_vector,
                                    mu_vector,
                                    eta_vector,
                                    correlation_structure,
                                    alpha_vector,
                                    phi);
  case MethodCode::pgee_jeffreys:
    return update_beta_jeffreys_cc(y_vector,
                                   model_matrix,
                                   id_vector,
                                   clusters,
                                   repeated_vector,
                                   weights_vector,
                                   link,
                                   family,
                                   beta_vector,
                                   mu_vector,
                                   eta_vector,
                                   correlation_structure,
                                   alpha_vector,
                                   phi,
                                   jeffreys_power);
  }
  Rcpp::stop("update_beta_cc: unsupported method.");
}
//==============================================================================


} // namespace

//============================ fitting function ================================
// [[Rcpp::export]]
Rcpp::List fit_geesolver_cc(const arma::vec& y_vector,
                            const arma::mat& model_matrix,
                            const arma::vec& id_vector,
                            const arma::vec& repeated_vector,
                            const arma::vec& weights_vector,
                            const char* link,
                            const char* family,
                            arma::vec beta_vector,
                            const arma::vec& offset,
                            const int maxiter,
                            const double tolerance,
                            const int step_maxiter,
                            const double step_multiplier,
                            const double jeffreys_power,
                            const char* method,
                            int use_params,
                            arma::vec alpha_vector,
                            const int alpha_fixed,
                            const char* correlation_structure,
                            const int mdependence,
                            double phi,
                            const int phi_fixed,
                            const int hold_nuisance) {
  const LinkCode link_code = parse_link(link);
  const FamilyCode family_code = parse_family(family);
  const arma::uword params_no = model_matrix.n_cols;
  use_params = (use_params != 0) ? static_cast<int>(params_no) : 0;
  const auto clusters = clusters_from_sorted_id(id_vector);
  // Iterates are stored in a growable container so that a large maxiter does
  // not allocate a p x (maxiter + 1) matrix up front.
  std::vector<arma::vec> beta_history;
  beta_history.reserve(static_cast<size_t>(std::min(maxiter, 100)) + 1);
  beta_history.push_back(beta_vector);
  arma::vec stepsize_vector(params_no, arma::fill::zeros);
  arma::vec criterion_vector(maxiter, arma::fill::zeros);
  arma::vec eta_vector = model_matrix * beta_vector + offset;
  if (!geer::is_valid_eta(link_code, eta_vector)) {
    Rcpp::stop(
      "invalid initial linear predictor: please try different starting values for beta."
    );
  }
  arma::vec mu_vector = geer::inverse_link(link_code, eta_vector);
  if (!geer::is_valid_mu(family_code, mu_vector)) {
    Rcpp::stop(
      "invalid initial fitted values: please try different starting values for beta."
    );
  }
  arma::vec pearson_residuals_vector =
    get_pearson_residuals(family_code,
                          y_vector,
                          mu_vector,
                          weights_vector);
  if (phi_fixed == 0) {
    phi = get_phi_hat(pearson_residuals_vector, use_params);
  }
  if (alpha_fixed == 0) {
    alpha_vector = get_alpha_hat(correlation_structure,
                                 pearson_residuals_vector,
                                 clusters,
                                 repeated_vector,
                                 phi,
                                 use_params,
                                 mdependence);
  }
  auto refresh_state = [&](const arma::vec& beta) {
    eta_vector = model_matrix * beta + offset;
    mu_vector = geer::inverse_link(link_code, eta_vector);
    pearson_residuals_vector =
      get_pearson_residuals(family_code,
                            y_vector,
                            mu_vector,
                            weights_vector);
    if (phi_fixed == 0 && hold_nuisance == 0) {
      phi = get_phi_hat(pearson_residuals_vector, use_params);
    }
    if (alpha_fixed == 0 && hold_nuisance == 0) {
      alpha_vector = get_alpha_hat(correlation_structure,
                                   pearson_residuals_vector,
                                   clusters,
                                   repeated_vector,
                                   phi,
                                   use_params,
                                   mdependence);
    }
  };
  std::string failure_message;
  arma::vec eta_trial_vector;
  arma::vec mu_trial_vector;
  for (int i = 1; i < maxiter + 1; ++i) {
    if (i == 1) {
      stepsize_vector =
        update_beta_cc(y_vector,
                       model_matrix,
                       id_vector,
                       clusters,
                       repeated_vector,
                       weights_vector,
                       link,
                       family,
                       beta_vector,
                       mu_vector,
                       eta_vector,
                       correlation_structure,
                       alpha_vector,
                       phi,
                       jeffreys_power,
                       method) - beta_vector;
      if (!stepsize_vector.is_finite()) {
        Rcpp::stop("the Newton step at the starting values contains "
                     "non-finite values: please try different starting "
                     "values for beta.");
      }
    }
    double criterion_inner = arma::norm(stepsize_vector, "inf");
    arma::vec beta_vector_inner = beta_vector;
    arma::vec stepsize_vector_inner = stepsize_vector;
    arma::vec beta_vector_new;
    arma::vec beta_vector_new_inner;
    bool valid_candidate_found = false;
    bool numerical_failure = false;
    int last_failure = 0;
    for (int j = 1; j < step_maxiter + 1; ++j) {
      beta_vector_new_inner =
        beta_vector_inner +
        step_multiplier * std::pow(0.5, j - 1) * stepsize_vector_inner;
      eta_trial_vector = model_matrix * beta_vector_new_inner + offset;
      if (!geer::is_valid_eta(link_code, eta_trial_vector)) {
        last_failure = 1;
        continue;
      }
      mu_trial_vector = geer::inverse_link(link_code, eta_trial_vector);
      if (!geer::is_valid_mu(family_code, mu_trial_vector)) {
        last_failure = 2;
        continue;
      }
      valid_candidate_found = true;
      eta_vector.swap(eta_trial_vector);
      mu_vector.swap(mu_trial_vector);
      try {
        pearson_residuals_vector =
          get_pearson_residuals(family_code,
                                y_vector,
                                mu_vector,
                                weights_vector);
        if (phi_fixed == 0 && hold_nuisance == 0) {
          phi = get_phi_hat(pearson_residuals_vector, use_params);
        }
        if (alpha_fixed == 0 && hold_nuisance == 0) {
          alpha_vector = get_alpha_hat(correlation_structure,
                                       pearson_residuals_vector,
                                       clusters,
                                       repeated_vector,
                                       phi,
                                       use_params,
                                       mdependence);
        }
        stepsize_vector_inner =
          update_beta_cc(y_vector,
                         model_matrix,
                         id_vector,
                         clusters,
                         repeated_vector,
                         weights_vector,
                         link,
                         family,
                         beta_vector_new_inner,
                         mu_vector,
                         eta_vector,
                         correlation_structure,
                         alpha_vector,
                         phi,
                         jeffreys_power,
                         method) - beta_vector_new_inner;
        if (!stepsize_vector_inner.is_finite()) {
          throw std::runtime_error(
            "the Newton step contains non-finite values");
        }
      } catch (const std::bad_alloc&) {
        throw;
      } catch (const std::exception& e) {
        numerical_failure = true;
        failure_message = e.what();
        break;
      }

      beta_vector_new = beta_vector_new_inner;
      beta_vector_inner = beta_vector_new_inner;
      const double criterion_candidate =
        arma::norm(stepsize_vector_inner, "inf");
      if (criterion_inner > criterion_candidate) {
        break;
      }
    }
    if (numerical_failure || !valid_candidate_found) {
      if (!numerical_failure) {
        failure_message = (last_failure == 2) ?
          "every step-halving attempt produced invalid fitted values" :
          "every step-halving attempt produced an invalid linear predictor";
      }
      if (valid_candidate_found) {
        refresh_state(beta_vector);
      }
      beta_history.push_back(beta_vector);
      criterion_vector(i - 1) = arma::datum::inf;
      break;
    }
    criterion_vector(i - 1) = arma::norm(stepsize_vector_inner, "inf");
    beta_vector = beta_vector_new;
    beta_history.push_back(beta_vector);
    stepsize_vector = stepsize_vector_inner;
    if (criterion_vector(i - 1) <= tolerance) {
      break;
    }
  }
  arma::mat beta_hat_matrix(params_no, beta_history.size());
  for (arma::uword col = 0; col < beta_history.size(); ++col) {
    beta_hat_matrix.col(col) = beta_history[col];
  }
  Rcpp::List cov_matrices =
    get_covariance_matrices_cc(y_vector,
                               model_matrix,
                               id_vector,
                               repeated_vector,
                               weights_vector,
                               link,
                               family,
                               mu_vector,
                               eta_vector,
                               correlation_structure,
                               alpha_vector,
                               phi);
  Rcpp::List result;
  result["beta_hat"] = beta_vector;
  result["beta_mat"] = beta_hat_matrix;
  result["alpha"] = alpha_vector;
  result["phi"] = phi;
  result["naive_covariance"] = cov_matrices["naive_covariance"];
  result["robust_covariance"] = cov_matrices["robust_covariance"];
  result["bc_covariance"] = cov_matrices["bc_covariance"];
  result["criterion"] = criterion_vector;
  result["eta"] = eta_vector;
  result["residuals"] = y_vector - mu_vector;
  result["fitted"] = mu_vector;
  result["offset"] = offset;
  result["failure"] = failure_message;
  return result;
}
//==============================================================================
