#include <RcppArmadillo.h>


//============================ estimate marginalized odds ratio structure ======
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector get_marginalized_odds_ratios(const arma::vec & response_vector,
                                                 const arma::vec & id_vector,
                                                 const arma::vec & repeated_vector,
                                                 const double & adding_constant,
                                                 Rcpp::String or_structure){
  int sample_size = max(id_vector);
  int cluster_size_max = max(repeated_vector);
  arma::mat wide_responses_matrix = arma::zeros(sample_size, cluster_size_max);
  for(int i = 1; i < sample_size + 1; i++) {
    arma::uvec id_vector_i = find(id_vector == i);
    arma::vec response_vector_i = response_vector(id_vector_i);
    arma::vec repeated_vector_i = repeated_vector(id_vector_i);
    int cluster_size_i = repeated_vector_i.n_elem;
    for(int j = 1; j < cluster_size_i + 1; j++) {
      wide_responses_matrix(i - 1, repeated_vector_i(j - 1) - 1) =
        response_vector_i(j - 1) + 1.0;
    }
  }
  double pairs_no = cluster_size_max * (cluster_size_max - 1.0) / 2.0;
  arma::vec counts = arma::zeros(4.0 * pairs_no);
  for(int l=1; l<sample_size+1; l++) {
    int k = 1;
    for(int categ_one = 1; categ_one < 3; categ_one++) {
      for(int categ_two = 1; categ_two < 3; categ_two++) {
        for(int i = 1; i < cluster_size_max; i++){
          if(wide_responses_matrix(l - 1, i - 1) == categ_one) {
            for(int j = i + 1; j < cluster_size_max + 1; j++) {
              if(wide_responses_matrix(l - 1, j - 1) == categ_two) {
                counts(k - 1) += 1;
              }
              k += 1;
            }
          } else {
            k += cluster_size_max - i;
          }
        }
      }
    }
  }
  counts += adding_constant;
  Rcpp::NumericVector ans(pairs_no);
  for(int i=1; i<pairs_no+1; i++) {
    ans(i-1) = (counts(i - 1) * counts(i - 1 + 3 * pairs_no)) /
      (counts(i - 1 + pairs_no) * counts(i - 1 + 2 * pairs_no));
  }
  Rcpp::String exchangeable("exchangeable");
  if (or_structure == exchangeable) {
    ans = log(ans);
    ans = Rcpp::rep(mean(ans), pairs_no);
    ans = exp(ans);
  }
  return ans;
}
//==============================================================================


//============================ subject specific odds ratios ====================
// [[Rcpp::export]]
arma::vec get_subject_specific_odds_ratios(arma::vec repeated_vector_i,
                                           int cluster_size_max,
                                           arma::vec odds_ratios_vector) {
  int cluster_size_i = repeated_vector_i.n_elem;
  double pairs_no_i = cluster_size_i * (cluster_size_i - 1)/2;
  arma::vec ans = arma::zeros(pairs_no_i);
  double pos = 2;
  int k = 1;
  for(int i = 1; i < cluster_size_i; i++){
    double row_index = repeated_vector_i(i - 1);
    for(int j = i + 1; j < cluster_size_i + 1; j++){
      pos = (row_index - 1) * (cluster_size_max - row_index / 2) +
        (repeated_vector_i(j - 1) - row_index) - 1;
      ans(k - 1) = odds_ratios_vector(pos);
      k = k + 1;
    }
  }
  return ans;
}
//==============================================================================


//============================ bivariate distribution ==========================
// [[Rcpp::export]]
double get_bivariate_distribution(const double & row_prob,
                                  const double & col_prob,
                                  const double & odds_ratio) {
  double ans = row_prob * col_prob;
  double tol = 0.000001;
  double odds_ratio_new = odds_ratio;
  if (row_prob > (1 - tol) || col_prob > (1.0 - tol)) {
    odds_ratio_new = 1;
  }
  if (row_prob < (tol) || col_prob < (tol)) {
    odds_ratio_new = 1;
  }
  if (odds_ratio_new == 1) {
    return ans;
  }
  double f_value = 1 - (1 - odds_ratio_new) * (row_prob + col_prob);
  double root_value =
    pow(f_value, 2) - 4 * odds_ratio_new * (odds_ratio_new - 1) * ans;
  ans = (f_value - sqrt(root_value)) / (2 * (odds_ratio_new - 1));
  return ans;
}
//==============================================================================


//============================ weight matrix ===================================
// [[Rcpp::export]]
arma::mat get_weight_matrix_or(const arma::vec & mu_vector,
                               const arma::vec & odds_ratios_vector,
                               const arma::vec & weights_vector) {
  // size of the probability vector
  int cluster_size = mu_vector.n_elem;
  // initializing the weight matrix V_i
  arma::mat ans = arma::diagmat(mu_vector % (1 - mu_vector));
  // calculating V_i if size > 1
  if(cluster_size > 1) {
    int k = 1;
    for(int i = 1; i < cluster_size; i++) {
      for(int j = i + 1; j < cluster_size + 1; j++) {
        ans(i - 1, j - 1) =
          get_bivariate_distribution(mu_vector(i - 1),
                                     mu_vector(j - 1),
                                     odds_ratios_vector(k - 1)) -
                                       mu_vector(i - 1) * mu_vector(j - 1);
        k += 1;
      }
    }
    ans = arma::symmatu(ans);
  }
  arma::vec weights_sq_inv_vector = 1/sqrt(weights_vector);
  ans = ans % (weights_sq_inv_vector * trans(weights_sq_inv_vector));
  return ans;
}
//==============================================================================


//============================ inverse weight matrix ===========================
// [[Rcpp::export]]
arma::mat get_weight_matrix_inverse_or(const arma::vec & mu_vector,
                                       const arma::vec & odds_ratios_vector,
                                       const arma::vec & weights_vector) {
  arma::mat weight_matrix = get_weight_matrix_or(mu_vector,
                                                 odds_ratios_vector,
                                                 weights_vector);
  arma::mat ans = arma::pinv(weight_matrix);
  return ans;
}
//==============================================================================


//============================ first derivative wrt row probability ============
// [[Rcpp::export]]
double get_bivariate_distribution_murow(const double & row_prob,
                                        const double & col_prob,
                                        const double & odds_ratio) {
  double f_value = 1 - (1 - odds_ratio) * (row_prob + col_prob);
  double biv_dis = get_bivariate_distribution(row_prob, col_prob, odds_ratio);
  double num = row_prob + col_prob - odds_ratio * (row_prob - col_prob) - 1;
  double den = f_value - 2 * (odds_ratio - 1) * biv_dis;
  double ans = 0.5 * (1 + num/den);
  return ans;
}
//==============================================================================


//============================ second derivative wrt row probability ===========
// [[Rcpp::export]]
double get_bivariate_distribution_murow2(const double & row_prob,
                                         const double & col_prob,
                                         const double & odds_ratio) {
  double f_value = 1 - (1 - odds_ratio) * (row_prob + col_prob);
  double biv_dis = get_bivariate_distribution(row_prob, col_prob, odds_ratio);
  double num = 2 * odds_ratio * (odds_ratio - 1) * col_prob * (col_prob - 1);
  double den = pow(f_value - 2 * (odds_ratio - 1) * biv_dis, 3);
  double ans = num/den;
  return ans;
}
//==============================================================================


//============================ second derivative wrt row-col probabilities =====
// [[Rcpp::export]]
double get_bivariate_distribution_murowcol(double row_prob,
                                           double col_prob,
                                           double odds_ratio) {
  double f_value = 1 - (1 - odds_ratio) * (row_prob + col_prob);
  double biv_dis = get_bivariate_distribution(row_prob, col_prob, odds_ratio);
  double num = (f_value - 2 * (odds_ratio - 1) * row_prob * col_prob) * odds_ratio;
  double den = pow(f_value - 2 * (odds_ratio - 1) * biv_dis, 3);
  double ans = num/den;
  return ans;
}
//==============================================================================


//============================ derivatives g_matrix ============================
// [[Rcpp::export]]
arma::mat get_g_matrix(const arma::vec & mu_vector,
                       const arma::vec & odds_ratios_vector) {
  int cluster_size = mu_vector.n_elem;
  arma::mat ans = arma::zeros(cluster_size, cluster_size);
  if(cluster_size > 1) {
    for(int i = 1; i < cluster_size; i++) {
      for(int j = i + 1; j < cluster_size + 1; j++) {
        double k = j - 0.5 * i * (i + 1) + cluster_size * (i - 1);
        ans(i - 1, j - 1) =
          get_bivariate_distribution_murow(mu_vector(i - 1),
                                           mu_vector(j - 1),
                                           odds_ratios_vector(k - 1));
      }
    }
    for(int i = 2; i < cluster_size + 1; i++) {
      for(int j = 1; j < i; j++) {
        double k = i - 0.5 * j * (j + 1) + cluster_size * (j - 1);
        ans(i - 1, j - 1) =
          get_bivariate_distribution_murow(mu_vector(i - 1),
                                           mu_vector(j - 1),
                                           odds_ratios_vector(k - 1));
      }
    }
  }
  return ans;
}
//==============================================================================


//============================ joint probability distribution derivative =======
// [[Rcpp::export]]
arma::mat get_bivariate_distribution_mu(const arma::vec & mu_vector,
                                        const arma::vec & odds_ratios_vector) {
  int cluster_size = mu_vector.n_elem;
  arma::mat ans = arma::zeros(pow(cluster_size, 2), cluster_size);
  if(cluster_size > 1) {
    for(int r = 1; r < cluster_size + 1; r++ ) {
      for(int i = 1; i < cluster_size + 1; i++ ) {
        if(i != r) {
          double k =
            std::max(i, r) - 0.5 * std::min(i, r) * (std::min(i, r) + 1) +
            cluster_size * (std::min(i, r) - 1);
          ans((r - 1) * cluster_size + i - 1, i - 1) =
            get_bivariate_distribution_murow(mu_vector(i - 1),
                                             mu_vector(r - 1),
                                             odds_ratios_vector(k - 1));
          ans((r - 1) * cluster_size + i - 1, r - 1) +=
            get_bivariate_distribution_murow(mu_vector(r - 1),
                                             mu_vector(i - 1),
                                             odds_ratios_vector(k - 1));
        }
      }
    }
  }
  return ans;
}
//==============================================================================


//============================ derivative of g matrix ==========================
// [[Rcpp::export]]
arma::mat get_g_matrix_mu(const arma::vec & mu_vector,
                          const arma::vec & odds_ratios_vector) {
  int cluster_size = mu_vector.n_elem;
  arma::mat ans = arma::zeros(pow(cluster_size, 2), cluster_size);
  // calculating V_i if size > 1
  if(cluster_size > 1) {
    for(int r = 1; r < cluster_size + 1; r++) {
      for(int i = 1; i < cluster_size + 1; i++) {
        if(i != r) {
          double k =
            std::max(i, r) - 0.5 * std::min(i, r) * (std::min(i, r) + 1) +
            cluster_size * (std::min(i, r) - 1);
          ans((r - 1) * cluster_size + i - 1, r - 1) =
            get_bivariate_distribution_murowcol(mu_vector(i - 1),
                                                mu_vector(r - 1),
                                                odds_ratios_vector(k - 1));
          ans((r - 1) * cluster_size + i - 1, i - 1) =
            get_bivariate_distribution_murow2(mu_vector(i - 1),
                                              mu_vector(r - 1),
                                              odds_ratios_vector(k - 1));
        }
      }
    }
  }
  return ans;
}
//==============================================================================


//============================ derivative of weight matrix =====================
// [[Rcpp::export]]
arma::mat get_weight_matrix_mu_or(const arma::vec & mu_vector,
                                  const arma::vec & odds_ratios_vector,
                                  const arma::vec & weights_vector) {
  int cluster_size = mu_vector.n_elem;
  arma::mat ans = arma::zeros(pow(cluster_size, 2), cluster_size);
  if(cluster_size > 1) {
    for(int r = 1; r < cluster_size + 1; r++) {
      for(int i = 1; i < cluster_size + 1; i++) {
        ans((r - 1) * cluster_size + i - 1, i - 1) = - mu_vector(r - 1);
        ans((r - 1) * cluster_size + i - 1, r - 1) -= mu_vector(i - 1);
        if(i != r) {
          double k =
            std::max(i, r) - 0.5 * std::min(i, r) * (std::min(i, r) + 1) +
            cluster_size * (std::min(i, r) - 1);
          ans((r - 1) * cluster_size + i - 1, i - 1) +=
            get_bivariate_distribution_murow(mu_vector(i - 1),
                                             mu_vector(r - 1),
                                             odds_ratios_vector(k - 1));
          ans((r - 1) * cluster_size + i - 1, r - 1) +=
            get_bivariate_distribution_murow(mu_vector(r - 1),
                                             mu_vector(i - 1),
                                             odds_ratios_vector(k - 1));
        }
      }
    }
    for(int j = 1; j < cluster_size + 1; j++) {
      ans((j - 1) * cluster_size + j - 1, j - 1) += 1;
    }
  }
  arma::mat weights_sq_inverse_matrix = arma::diagmat(1/sqrt(weights_vector));
  ans = kron(weights_sq_inverse_matrix, weights_sq_inverse_matrix) * ans;
  return ans;
}
//==============================================================================
