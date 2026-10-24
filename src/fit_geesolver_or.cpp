#include "cluster_utils.h"
#include "covariance_matrices.h"
#include "link_functions.h"
#include "method_codes.h"
#include "utils.h"
#include "working_covariance_or.h"
#include <algorithm>
#include <cmath>
#include <new>
#include <stdexcept>
#include <string>
#include <vector>


namespace {

//============================ update beta - gee OR ============================
arma::vec update_beta_gee_or(const arma::vec& y_vector,
                             const arma::mat& model_matrix,
                             const arma::vec& id_vector,
                             const std::vector<Cluster>& clusters,
                             const arma::vec& repeated_vector,
                             const arma::vec& weights_vector,
                             const char* link,
                             const arma::vec& beta_vector,
                             const arma::vec& mu_vector,
                             const arma::vec& eta_vector,
                             const arma::vec& alpha_vector) {
  const LinkCode link_code = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max = static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = link_derivative_1(link_code, eta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  arma::mat v_matrix_inverse_d_matrix_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.first_row;
      const arma::uword last_row = cluster.after_last_row - 1;
      const ClusterWorkspace workspace = build_cluster_data_or(
        cluster, model_matrix, mu_vector, repeated_vector, weights_vector,
        delta_vector, s_vector, repeated_max, alpha_vector);
      const ClusterData& cluster_data = workspace.data;
      const arma::mat& d_matrix_i = cluster_data.d_matrix;
      const arma::mat& v_matrix_i = cluster_data.v_matrix;
      v_matrix_inverse_d_matrix_i =
        solve_chol_or_lu_mat(v_matrix_i, d_matrix_i);
      d_matrix_trans_v_matrix_inverse_i = v_matrix_inverse_d_matrix_i.t();
      information_matrix += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector.subvec(first_row, last_row);
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("update_beta_gee_or", cluster_index,
                                   id_vector[cluster.first_row], e);
    }
  }
  symmetrize_if_close(information_matrix, "update_beta_gee_or");
  return beta_vector +
         solve_chol_or_lu_vec(information_matrix, u_vector,
                              "update_beta_gee_or: naive information matrix");
}
//==============================================================================


//============================ update beta - naive OR ==========================
arma::vec update_beta_naive_or(const arma::vec& y_vector,
                               const arma::mat& model_matrix,
                               const arma::vec& id_vector,
                               const std::vector<Cluster>& clusters,
                               const arma::vec& repeated_vector,
                               const arma::vec& weights_vector,
                               const char* link,
                               const arma::vec& beta_vector,
                               const arma::vec& mu_vector,
                               const arma::vec& eta_vector,
                               const arma::vec& alpha_vector) {
  const LinkCode link_code = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max = static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat lambda_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = link_derivative_1(link_code, eta_vector);
  const arma::vec delta_star_vector =
    link_derivative_2(link_code, eta_vector) / arma::square(delta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  arma::mat d_matrix_trans_i;
  arma::mat v_matrix_inverse_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::mat weights_matrix_sq_inverse_i;
  arma::mat v_matrix_tilde_inverse_i;
  arma::mat g_matrix_i;
  arma::mat identity_matrix_i;
  arma::mat vecdiag_matrix_delta_star_matrix_i;
  arma::vec v_matrix_tilde_inverse_mu_vector_i;
  arma::mat h_epsilon_trans_i;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.first_row;
      const arma::uword last_row = cluster.after_last_row - 1;
      const arma::uword cluster_size =
        cluster.after_last_row - cluster.first_row;
      const arma::vec mu_vector_i = mu_vector.subvec(first_row, last_row);
      const arma::vec weights_vector_i = weights_vector.subvec(first_row, last_row);
      const arma::vec s_vector_i = s_vector.subvec(first_row, last_row);
      const ClusterWorkspace workspace = build_cluster_data_or(
        cluster, model_matrix, mu_vector, repeated_vector, weights_vector,
        delta_vector, s_vector, repeated_max, alpha_vector);
      const ClusterData& cluster_data = workspace.data;
      const arma::mat& d_matrix_i = cluster_data.d_matrix;
      const arma::mat& v_matrix_i = cluster_data.v_matrix;
      const arma::vec& odds_ratios_vector_i = workspace.odds_ratios;
      d_matrix_trans_i = d_matrix_i.t();
      v_matrix_inverse_i = solve_chol_or_lu_mat(
        v_matrix_i, arma::eye(cluster_size, cluster_size));
      d_matrix_trans_v_matrix_inverse_i = d_matrix_trans_i * v_matrix_inverse_i;
      information_matrix += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      const arma::vec u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector_i;
      u_vector += u_vector_i;
      weights_matrix_sq_inverse_i = arma::diagmat(1.0 / arma::sqrt(weights_vector_i));
      v_matrix_tilde_inverse_i = v_matrix_inverse_i * weights_matrix_sq_inverse_i;
      g_matrix_i = get_g_matrix_or(mu_vector_i, odds_ratios_vector_i);
      identity_matrix_i = arma::eye(cluster_size, cluster_size);
      vecdiag_matrix_delta_star_matrix_i =
        vecdiag_right_diag(delta_star_vector.subvec(first_row, last_row));
      v_matrix_tilde_inverse_mu_vector_i =
        v_matrix_tilde_inverse_i * mu_vector_i;
      h_epsilon_trans_i =
        vecdiag_matrix_delta_star_matrix_i +
        (arma::vectorise(v_matrix_tilde_inverse_i.t()) * mu_vector_i.t() -
        kronecker_left_identity_vecdiag(v_matrix_tilde_inverse_i) * g_matrix_i -
        kronecker_left_identity_vecdiag(
          v_matrix_tilde_inverse_i * (g_matrix_i.t() + identity_matrix_i)
        ) +
          arma::kron(v_matrix_tilde_inverse_mu_vector_i, identity_matrix_i)) *
          weights_matrix_sq_inverse_i;
      add_kronecker_self_t_s_d(
        lambda_matrix,
        d_matrix_i,
        arma::kron(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i.t()) -
          kronecker_sum_same(v_matrix_inverse_i) * vecdiag_matrix_delta_star_matrix_i -
          arma::kron(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i) +
          h_epsilon_trans_i * v_matrix_inverse_i +
          arma::kron(v_matrix_inverse_i, v_matrix_inverse_i) *
            get_v_matrix_mu_or(mu_vector_i,
                               odds_ratios_vector_i,
                               weights_vector_i),
        -1.0);
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("update_beta_naive_or", cluster_index,
                                   id_vector[cluster.first_row], e);
    }
  }
  symmetrize_if_close(information_matrix, "update_beta_naive_or");
  arma::vec lambda_vector(params_no, arma::fill::zeros);
  for (arma::uword r = 0; r < params_no; ++r) {
    const arma::mat block =
      lambda_matrix.rows(r * params_no, (r + 1) * params_no - 1);
    const arma::mat solved_block =
      solve_chol_or_lu_mat(information_matrix, block,
                           "update_beta_naive_or: naive information matrix");
    lambda_vector[r] = 0.5 * arma::trace(solved_block);
  }
  const arma::vec update_step =
    solve_chol_or_lu_vec(information_matrix, u_vector + lambda_vector,
                         "update_beta_naive_or: naive information matrix");
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - robust OR =========================
arma::vec update_beta_robust_or(const arma::vec& y_vector,
                                const arma::mat& model_matrix,
                                const arma::vec& id_vector,
                                const std::vector<Cluster>& clusters,
                                const arma::vec& repeated_vector,
                                const arma::vec& weights_vector,
                                const char* link,
                                const arma::vec& beta_vector,
                                const arma::vec& mu_vector,
                                const arma::vec& eta_vector,
                                const arma::vec& alpha_vector) {
  const LinkCode link_code = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max = static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat partial_derivatives_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat second_derivatives_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = link_derivative_1(link_code, eta_vector);
  const arma::vec delta_star_vector =
    link_derivative_2(link_code, eta_vector) / arma::square(delta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  arma::mat d_matrix_trans_i;
  arma::mat v_matrix_inverse_i;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::mat weights_matrix_sq_inverse_i;
  arma::mat v_matrix_tilde_inverse_i;
  arma::mat g_matrix_i;
  arma::mat identity_matrix_i;
  arma::mat vecdiag_matrix_delta_star_matrix_i;
  arma::vec v_matrix_tilde_inverse_mu_vector_i;
  arma::mat h_epsilon_trans_i;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.first_row;
      const arma::uword last_row = cluster.after_last_row - 1;
      const arma::uword cluster_size =
        cluster.after_last_row - cluster.first_row;
      const arma::vec mu_vector_i = mu_vector.subvec(first_row, last_row);
      const arma::vec weights_vector_i = weights_vector.subvec(first_row, last_row);
      const arma::vec s_vector_i = s_vector.subvec(first_row, last_row);
      const ClusterWorkspace workspace = build_cluster_data_or(
        cluster, model_matrix, mu_vector, repeated_vector, weights_vector,
        delta_vector, s_vector, repeated_max, alpha_vector);
      const ClusterData& cluster_data = workspace.data;
      const arma::mat& d_matrix_i = cluster_data.d_matrix;
      const arma::mat& v_matrix_i = cluster_data.v_matrix;
      const arma::vec& odds_ratios_vector_i = workspace.odds_ratios;
      d_matrix_trans_i = d_matrix_i.t();
      v_matrix_inverse_i = solve_chol_or_lu_mat(
        v_matrix_i, arma::eye(cluster_size, cluster_size));
      d_matrix_trans_v_matrix_inverse_i = d_matrix_trans_i * v_matrix_inverse_i;
      information_matrix += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      const arma::vec u_vector_i = d_matrix_trans_v_matrix_inverse_i * s_vector_i;
      u_vector += u_vector_i;
      meat_matrix += u_vector_i * u_vector_i.t();
      weights_matrix_sq_inverse_i = arma::diagmat(1.0 / arma::sqrt(weights_vector_i));
      v_matrix_tilde_inverse_i = v_matrix_inverse_i * weights_matrix_sq_inverse_i;
      g_matrix_i = get_g_matrix_or(mu_vector_i, odds_ratios_vector_i);
      identity_matrix_i = arma::eye(cluster_size, cluster_size);
      vecdiag_matrix_delta_star_matrix_i =
        vecdiag_right_diag(delta_star_vector.subvec(first_row, last_row));
      v_matrix_tilde_inverse_mu_vector_i =
        v_matrix_tilde_inverse_i * mu_vector_i;
      h_epsilon_trans_i =
        vecdiag_matrix_delta_star_matrix_i +
        (arma::vectorise(v_matrix_tilde_inverse_i.t()) * mu_vector_i.t() -
        kronecker_left_identity_vecdiag(v_matrix_tilde_inverse_i) * g_matrix_i -
        kronecker_left_identity_vecdiag(
          v_matrix_tilde_inverse_i * (g_matrix_i.t() + identity_matrix_i)
        ) +
          arma::kron(v_matrix_tilde_inverse_mu_vector_i, identity_matrix_i)) *
          weights_matrix_sq_inverse_i;
      partial_derivatives_matrix +=
        kronecker_self_t_vec(
          d_matrix_i,
          h_epsilon_trans_i * v_matrix_inverse_i * s_vector_i) *
        u_vector_i.t();
      add_kronecker_self_t_s_d(
        second_derivatives_matrix,
        d_matrix_i,
        arma::kron(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i.t()) -
          kronecker_sum_same(v_matrix_inverse_i) * vecdiag_matrix_delta_star_matrix_i -
          arma::kron(v_matrix_tilde_inverse_mu_vector_i, v_matrix_tilde_inverse_i) -
          h_epsilon_trans_i * v_matrix_inverse_i +
          arma::kron(v_matrix_inverse_i, v_matrix_inverse_i) *
            get_v_matrix_mu_or(mu_vector_i,
                               odds_ratios_vector_i,
                               weights_vector_i));
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("update_beta_robust_or", cluster_index,
                                   id_vector[cluster.first_row], e);
    }
  }
  symmetrize_if_close(information_matrix, "update_beta_robust_or");
  const arma::mat naive_covariance_meat_matrix =
    solve_chol_or_lu_mat(information_matrix, meat_matrix,
                         "update_beta_robust_or: naive information matrix");
  const arma::mat robust_matrix =
    solve_chol_or_lu_mat(information_matrix, naive_covariance_meat_matrix.t(),
                         "update_beta_robust_or: naive information matrix");
  arma::vec lambda_vector(params_no, arma::fill::zeros);
  for (arma::uword r = 0; r < params_no; ++r) {
    const arma::mat first_block =
      partial_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
    const arma::mat second_block =
      second_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
    const arma::mat solved_first_block =
      solve_chol_or_lu_mat(information_matrix, first_block,
                           "update_beta_robust_or: naive information matrix");
    lambda_vector[r] =
      -(arma::trace(solved_first_block) +
        0.5 * arma::trace(robust_matrix * second_block));
  }
  const arma::vec update_step =
    solve_chol_or_lu_vec(information_matrix, u_vector + lambda_vector,
                         "update_beta_robust_or: naive information matrix");
  return beta_vector + update_step;
}
//==============================================================================


//============================ derivative of epsilon - OR ======================
arma::mat get_epsilon_derivative_or(const arma::vec& mu_vector,
                                    const arma::vec& s_vector,
                                    const arma::mat& v_matrix_inverse,
                                    const arma::vec& weights_sq_inverse_vector,
                                    const arma::mat& g_matrix,
                                    const arma::mat& g_matrix_mu,
                                    const arma::vec& delta_star_vector,
                                    const arma::vec& delta_tilde_star_vector,
                                    const arma::mat& epsilon_matrix) {
  const arma::uword m = mu_vector.n_elem;
  arma::mat v_tilde = v_matrix_inverse;
  v_tilde.each_row() %= weights_sq_inverse_vector.t();
  const arma::vec q_vector = v_tilde.t() * s_vector;
  const arma::vec z_vector = v_matrix_inverse * s_vector;
  arma::mat h_matrix(m, m);
  for (arma::uword b = 0; b < m; ++b) {
    h_matrix.col(b) = g_matrix.row(b).t() - mu_vector;
    h_matrix(b, b) += 0.5;
  }
  const arma::mat r_matrix = v_tilde * h_matrix;
  const arma::vec g_vector = h_matrix.t() * q_vector;
  auto dp = [&](arma::uword a, arma::uword b, arma::uword c) {
    return -(v_tilde(a, c) * r_matrix(b, c) + r_matrix(a, c) * v_tilde(b, c));
  };
  arma::mat dq_matrix(m, m);
  for (arma::uword b = 0; b < m; ++b) {
    for (arma::uword c = 0; c < m; ++c) {
      dq_matrix(b, c) =
        weights_sq_inverse_vector[b] *
          (-q_vector[c] * r_matrix(b, c) - g_vector[c] * v_tilde(b, c)) -
        v_tilde(c, b);
    }
  }
  arma::mat dg_matrix(m, m);
  for (arma::uword b = 0; b < m; ++b) {
    for (arma::uword c = 0; c < m; ++c) {
      double value = arma::dot(h_matrix.col(b), dq_matrix.col(c));
      for (arma::uword j = 0; j < m; ++j) {
        value += (g_matrix_mu(b + j * m, c) - (j == c ? 1.0 : 0.0)) * q_vector[j];
      }
      dg_matrix(b, c) = value;
    }
  }
  arma::mat result(m * m, m, arma::fill::zeros);
  for (arma::uword a = 0; a < m; ++a) {
    for (arma::uword b = 0; b < m; ++b) {
      for (arma::uword c = 0; c < m; ++c) {
        const double dp_abc = dp(a, b, c);
        double value = -dp_abc;
        if (a == b) {
          const double dz = epsilon_matrix(a, c) -
                            (a == c ? delta_star_vector[a] * z_vector[a] : 0.0);
          value += (a == c ? delta_tilde_star_vector[a] * z_vector[a] : 0.0) +
                   delta_star_vector[a] * dz;
        }
        value -= dp_abc * weights_sq_inverse_vector[b] * g_vector[b] +
                 v_tilde(a, b) * dg_matrix(b, c);
        double d_vh = 0.0;
        for (arma::uword j = 0; j < m; ++j) {
          d_vh += dp(a, j, c) * weights_sq_inverse_vector[j] * h_matrix(j, b) +
                  v_tilde(a, j) * (g_matrix_mu(b + j * m, c) - (j == c ? 1.0 : 0.0));
        }
        value -= d_vh * q_vector[b] + r_matrix(a, b) * dq_matrix(b, c);
        result(a * m + b, c) = value;
      }
    }
  }
  return result;
}
//==============================================================================


//============================ update beta - empirical OR ======================
arma::vec update_beta_empirical_or(const arma::vec& y_vector,
                                   const arma::mat& model_matrix,
                                   const arma::vec& id_vector,
                                   const std::vector<Cluster>& clusters,
                                   const arma::vec& repeated_vector,
                                   const arma::vec& weights_vector,
                                   const char* link,
                                   const arma::vec& beta_vector,
                                   const arma::vec& mu_vector,
                                   const arma::vec& eta_vector,
                                   const arma::vec& alpha_vector) {
  const LinkCode link_code = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max = static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat partial_derivatives_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat second_derivatives_matrix(params_no * params_no, params_no, arma::fill::zeros);
  arma::mat observed_fisher_info_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat meat_matrix(params_no, params_no, arma::fill::zeros);
  const arma::vec delta_vector = link_derivative_1(link_code, eta_vector);
  const arma::vec mueta2_vector = link_derivative_2(link_code, eta_vector);
  const arma::vec delta_star_vector =
    mueta2_vector / arma::square(delta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  const arma::vec delta_tilde_star_vector =
    (delta_vector % link_derivative_3(link_code, eta_vector) - 2.0 * arma::square(mueta2_vector)) /
    arma::pow(delta_vector, 4.0);
  arma::mat d_matrix_trans_i;
  arma::mat v_matrix_inverse_i;
  arma::mat delta_star_matrix_i;
  arma::mat w_matrix_i;
  arma::mat identity_matrix_i;
  arma::mat v_matrix_inverse_derivative_i;
  arma::mat epsilon_matrix_i;
  arma::mat observed_fisher_info_matrix_i;
  arma::mat epsilon_matrix_transpose_derivative_i;
  arma::mat second_derivatives_matrix_terms12_i;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.first_row;
      const arma::uword last_row = cluster.after_last_row - 1;
      const arma::uword cluster_size =
        cluster.after_last_row - cluster.first_row;
      const arma::vec mu_vector_i = mu_vector.subvec(first_row, last_row);
      const arma::vec s_vector_i = s_vector.subvec(first_row, last_row);
      const arma::vec weights_vector_i = weights_vector.subvec(first_row, last_row);
      const ClusterWorkspace workspace = build_cluster_data_or(
        cluster, model_matrix, mu_vector, repeated_vector, weights_vector,
        delta_vector, s_vector, repeated_max, alpha_vector);
      const ClusterData& cluster_data = workspace.data;
      const arma::mat& d_matrix_i = cluster_data.d_matrix;
      const arma::mat& v_matrix_i = cluster_data.v_matrix;
      const arma::vec& odds_ratios_vector_i = workspace.odds_ratios;
      d_matrix_trans_i = d_matrix_i.t();
      v_matrix_inverse_i = solve_chol_or_lu_mat(
        v_matrix_i, arma::eye(cluster_size, cluster_size));
      const arma::vec u_vector_i = d_matrix_trans_i * v_matrix_inverse_i * s_vector_i;
      u_vector += u_vector_i;
      meat_matrix += u_vector_i * u_vector_i.t();
      delta_star_matrix_i = arma::diagmat(delta_star_vector.subvec(first_row, last_row));
      w_matrix_i = arma::diagmat(v_matrix_inverse_i * s_vector_i);
      identity_matrix_i = arma::eye(cluster_size, cluster_size);
      v_matrix_inverse_derivative_i =
        -arma::kron(v_matrix_inverse_i, v_matrix_inverse_i) *
        get_v_matrix_mu_or(mu_vector_i,
                           odds_ratios_vector_i,
                           weights_vector_i);
      epsilon_matrix_i =
        delta_star_matrix_i * w_matrix_i -
        v_matrix_inverse_i +
        arma::kron(s_vector_i.t(), identity_matrix_i) * v_matrix_inverse_derivative_i;
      observed_fisher_info_matrix_i =
        d_matrix_trans_i * epsilon_matrix_i * d_matrix_i;
      observed_fisher_info_matrix -= observed_fisher_info_matrix_i;
      partial_derivatives_matrix +=
        arma::vectorise(observed_fisher_info_matrix_i.t()) * u_vector_i.t();
      epsilon_matrix_transpose_derivative_i =
        get_epsilon_derivative_or(mu_vector_i,
                                  s_vector_i,
                                  v_matrix_inverse_i,
                                  1.0 / arma::sqrt(weights_vector_i),
                                  get_g_matrix_or(mu_vector_i, odds_ratios_vector_i),
                                  get_g_matrix_mu_or(mu_vector_i, odds_ratios_vector_i),
                                  delta_star_vector.subvec(first_row, last_row),
                                  delta_tilde_star_vector.subvec(first_row, last_row),
                                  epsilon_matrix_i);
      second_derivatives_matrix_terms12_i =
        kronecker_left_identity_vecdiag(epsilon_matrix_i * delta_star_matrix_i) +
        kronecker_identity_right_vecdiag(epsilon_matrix_i.t() * delta_star_matrix_i);
      add_kronecker_self_t_s_d(second_derivatives_matrix,
                          d_matrix_i,
                          second_derivatives_matrix_terms12_i +
                            epsilon_matrix_transpose_derivative_i);
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("update_beta_empirical_or", cluster_index,
                                   id_vector[cluster.first_row], e);
    }
  }
  arma::mat robust_matrix(params_no, params_no, arma::fill::zeros);
  arma::vec lambda_vector(params_no, arma::fill::zeros);
  const CholOrLuFactor fisher_factor(
    observed_fisher_info_matrix,
    "update_beta_empirical_or: observed Fisher information",
    false);
  const arma::mat fisher_inv_meat = fisher_factor.solve(meat_matrix);
  robust_matrix = fisher_factor.solve(arma::mat(fisher_inv_meat.t()));
  for (arma::uword r = 0; r < params_no; ++r) {
    const arma::mat first_block =
      partial_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
    const arma::mat second_block =
      second_derivatives_matrix.rows(r * params_no, (r + 1) * params_no - 1);
    const arma::mat solved_first_block =
      fisher_factor.solve(first_block);
    lambda_vector[r] =
      -(arma::trace(solved_first_block) +
        0.5 * arma::trace(robust_matrix * second_block));
  }
  const arma::vec update_step =
    fisher_factor.solve(arma::vec(u_vector + lambda_vector));
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - jeffreys OR =======================
arma::vec update_beta_jeffreys_or(const arma::vec& y_vector,
                                  const arma::mat& model_matrix,
                                  const arma::vec& id_vector,
                                  const std::vector<Cluster>& clusters,
                                  const arma::vec& repeated_vector,
                                  const arma::vec& weights_vector,
                                  const char* link,
                                  const arma::vec& beta_vector,
                                  const arma::vec& mu_vector,
                                  const arma::vec& eta_vector,
                                  const arma::vec& alpha_vector,
                                  const double jeffreys_power) {
  const LinkCode link_code = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const arma::uword repeated_max = static_cast<arma::uword>(arma::max(repeated_vector));
  arma::vec u_vector(params_no, arma::fill::zeros);
  arma::mat information_matrix(params_no, params_no, arma::fill::zeros);
  arma::mat naive_matrix_inverse_derivative(params_no * params_no,
                                            params_no,
                                            arma::fill::zeros);
  const arma::vec delta_vector = link_derivative_1(link_code, eta_vector);
  const arma::vec delta_star_vector =
    link_derivative_2(link_code, eta_vector) / arma::square(delta_vector);
  const arma::vec s_vector = y_vector - mu_vector;
  arma::mat d_matrix_trans_v_matrix_inverse_i;
  arma::mat v_matrix_inverse_i;
  for (arma::uword cluster_index = 0;
       cluster_index < clusters.size();
       ++cluster_index) {
    const auto& cluster = clusters[cluster_index];
    try {
      const arma::uword first_row = cluster.first_row;
      const arma::uword last_row = cluster.after_last_row - 1;
      const arma::uword cluster_size =
        cluster.after_last_row - cluster.first_row;
      const arma::vec mu_vector_i = mu_vector.subvec(first_row, last_row);
      const arma::vec weights_vector_i = weights_vector.subvec(first_row, last_row);
      const ClusterWorkspace workspace = build_cluster_data_or(
        cluster, model_matrix, mu_vector, repeated_vector, weights_vector,
        delta_vector, s_vector, repeated_max, alpha_vector);
      const ClusterData& cluster_data = workspace.data;
      const arma::mat& d_matrix_i = cluster_data.d_matrix;
      const arma::mat& v_matrix_i = cluster_data.v_matrix;
      const arma::vec& odds_ratios_vector_i = workspace.odds_ratios;
      v_matrix_inverse_i = solve_chol_or_lu_mat(
        v_matrix_i, arma::eye(cluster_size, cluster_size));
      d_matrix_trans_v_matrix_inverse_i = d_matrix_i.t() * v_matrix_inverse_i;
      information_matrix += d_matrix_trans_v_matrix_inverse_i * d_matrix_i;
      u_vector += d_matrix_trans_v_matrix_inverse_i * s_vector.subvec(first_row, last_row);
      add_kronecker_self_t_s_d(
        naive_matrix_inverse_derivative,
        d_matrix_i,
        kronecker_sum_same(
          v_matrix_inverse_i * arma::diagmat(delta_star_vector.subvec(first_row, last_row))
        ) *
          vecdiag_matrix(cluster_size) -
          arma::kron(v_matrix_inverse_i, v_matrix_inverse_i) *
            get_v_matrix_mu_or(mu_vector_i,
                               odds_ratios_vector_i,
                               weights_vector_i));
    } catch (const std::exception& e) {
      rethrow_with_cluster_context("update_beta_jeffreys_or", cluster_index,
                                   id_vector[cluster.first_row], e);
    }
  }
  symmetrize_if_close(information_matrix, "update_beta_jeffreys_or");
  const arma::mat naive_covariance =
    solve_chol_or_lu_mat(information_matrix, arma::eye(params_no, params_no),
                         "update_beta_jeffreys_or: naive information matrix");
  const arma::vec lambda_vector =
    jeffreys_power *
    naive_matrix_inverse_derivative.t() *
    arma::vectorise(naive_covariance);
  const arma::vec update_step =
    solve_chol_or_lu_vec(information_matrix, u_vector + lambda_vector,
                         "update_beta_jeffreys_or: naive information matrix");
  return beta_vector + update_step;
}
//==============================================================================


//============================ update beta - or ================================
arma::vec update_beta_or(const arma::vec& y_vector,
                         const arma::mat& model_matrix,
                         const arma::vec& id_vector,
                         const std::vector<Cluster>& clusters,
                         const arma::vec& repeated_vector,
                         const arma::vec& weights_vector,
                         const char* link,
                         const arma::vec& beta_vector,
                         const arma::vec& mu_vector,
                         const arma::vec& eta_vector,
                         const arma::vec& alpha_vector,
                         const double jeffreys_power,
                         const char* method) {
  switch (parse_method(method)) {
  case MethodCode::gee:
    return update_beta_gee_or(y_vector,
                              model_matrix,
                              id_vector,
                              clusters,
                              repeated_vector,
                              weights_vector,
                              link,
                              beta_vector,
                              mu_vector,
                              eta_vector,
                              alpha_vector);
  case MethodCode::brgee_naive:
    return update_beta_naive_or(y_vector,
                                model_matrix,
                                id_vector,
                                clusters,
                                repeated_vector,
                                weights_vector,
                                link,
                                beta_vector,
                                mu_vector,
                                eta_vector,
                                alpha_vector);
  case MethodCode::brgee_robust:
    return update_beta_robust_or(y_vector,
                                 model_matrix,
                                 id_vector,
                                 clusters,
                                 repeated_vector,
                                 weights_vector,
                                 link,
                                 beta_vector,
                                 mu_vector,
                                 eta_vector,
                                 alpha_vector);
  case MethodCode::brgee_empirical:
    return update_beta_empirical_or(y_vector,
                                    model_matrix,
                                    id_vector,
                                    clusters,
                                    repeated_vector,
                                    weights_vector,
                                    link,
                                    beta_vector,
                                    mu_vector,
                                    eta_vector,
                                    alpha_vector);
  case MethodCode::pgee_jeffreys:
    return update_beta_jeffreys_or(y_vector,
                                   model_matrix,
                                   id_vector,
                                   clusters,
                                   repeated_vector,
                                   weights_vector,
                                   link,
                                   beta_vector,
                                   mu_vector,
                                   eta_vector,
                                   alpha_vector,
                                   jeffreys_power);
  }
  Rcpp::stop("update_beta_or: unsupported method");
}
//==============================================================================


}  // namespace

//============================ fitting function ================================
// [[Rcpp::export]]
Rcpp::List fit_geesolver_or(const arma::vec& y_vector,
                            const arma::mat& model_matrix,
                            const arma::vec& id_vector,
                            const arma::vec& repeated_vector,
                            const arma::vec& weights_vector,
                            const char* link,
                            arma::vec beta_vector,
                            const arma::vec& offset,
                            const int maxiter,
                            const double tolerance,
                            const int step_maxiter,
                            const double step_multiplier,
                            const double jeffreys_power,
                            const char* method,
                            const arma::vec& alpha_vector) {
  reset_lu_fallback_warning();
  if (maxiter < 0 || step_maxiter < 0) {
    Rcpp::stop("'maxiter' and 'step_maxiter' must be nonnegative integers");
  }
  if (!model_matrix.is_finite() || !offset.is_finite() ||
      !weights_vector.is_finite() || !beta_vector.is_finite()) {
    Rcpp::stop("the model matrix, offset, weights and starting values must be finite");
  }
  const LinkCode link_code = parse_link(link);
  const arma::uword params_no = model_matrix.n_cols;
  const auto clusters = clusters_from_sorted_id(id_vector);
  for (arma::uword i = 0; i < y_vector.n_elem; ++i) {
    const double y_i = y_vector[i];
    if (!(y_i >= 0.0 && y_i <= 1.0)) {
      Rcpp::stop("the response must lie in [0, 1]; observation %d has value %g",
                 static_cast<int>(i + 1), y_i);
    }
  }
  std::vector<arma::vec> beta_history;
  beta_history.reserve(static_cast<size_t>(std::min(maxiter, 100)) + 1);
  beta_history.push_back(beta_vector);
  arma::vec stepsize_vector(params_no, arma::fill::zeros);
  arma::vec criterion_vector(maxiter, arma::fill::zeros);
  arma::vec eta_vector = model_matrix * beta_vector + offset;
  if (!is_valid_eta(link_code, eta_vector)) {
    Rcpp::stop(
      "the initial linear predictor is invalid; supply different 'beta_start' values"
    );
  }
  arma::vec mu_vector = inverse_link(link_code, eta_vector);
  if (!is_valid_mu(FamilyCode::binomial, mu_vector)) {
    Rcpp::stop(
      "the initial fitted values are invalid; supply different 'beta_start' values"
    );
  }
  auto refresh_state = [&](const arma::vec& beta) {
    eta_vector = model_matrix * beta + offset;
    mu_vector = inverse_link(link_code, eta_vector);
  };
  std::string failure_message;
  arma::vec eta_trial_vector;
  arma::vec mu_trial_vector;
  for (int i = 1; i < maxiter + 1; ++i) {
    if (i == 1) {
      stepsize_vector =
        update_beta_or(y_vector,
                       model_matrix,
                       id_vector,
                       clusters,
                       repeated_vector,
                       weights_vector,
                       link,
                       beta_vector,
                       mu_vector,
                       eta_vector,
                       alpha_vector,
                       jeffreys_power,
                       method) - beta_vector;
      if (!stepsize_vector.is_finite()) {
        Rcpp::stop("the Newton step at the starting values contains non-finite values; supply "
                   "different 'beta_start' values");
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
      if (!is_valid_eta(link_code, eta_trial_vector)) {
        last_failure = 1;
        continue;
      }
      mu_trial_vector = inverse_link(link_code, eta_trial_vector);
      if (!is_valid_mu(FamilyCode::binomial, mu_trial_vector)) {
        last_failure = 2;
        continue;
      }
      valid_candidate_found = true;
      eta_vector.swap(eta_trial_vector);
      mu_vector.swap(mu_trial_vector);
      try {
        stepsize_vector_inner =
          update_beta_or(y_vector,
                         model_matrix,
                         id_vector,
                         clusters,
                         repeated_vector,
                         weights_vector,
                         link,
                         beta_vector_new_inner,
                         mu_vector,
                         eta_vector,
                         alpha_vector,
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
    get_covariance_matrices_or(y_vector,
                               model_matrix,
                               id_vector,
                               repeated_vector,
                               weights_vector,
                               link,
                               mu_vector,
                               eta_vector,
                               alpha_vector);
  Rcpp::List result;
  result["beta_hat"] = beta_vector;
  result["beta_mat"] = beta_hat_matrix;
  result["alpha"] = alpha_vector;
  result["phi"] = 1.0;
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
