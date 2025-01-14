#ifndef UTILS_H
#define UTILS_H


arma::mat subset_matrix(arma::mat x, arma::vec y);

arma::mat kappa_matrix(int dimension);

arma::mat kronecker_sum_same(arma::mat x);

arma::mat kronecker_kappa1(arma::mat x);

arma::mat kronecker_kappa2(arma::mat x);

arma::mat kronecker_kappa3(arma::mat x);

arma::mat kronecker_left_identity_kappa(arma::mat x);

arma::mat kronecker_identity_right_kappa(arma::mat x);

arma::mat kappa_right(arma::mat x);

arma::mat kronecker_vector_identity(arma::vec x);

arma::mat kronecker_vector_identity(arma::vec x);

arma::mat kronecker_vector_matrix(arma::vec x, arma::mat y);


#endif
