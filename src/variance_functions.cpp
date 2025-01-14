#include <RcppArmadillo.h>
#include "utils.h"
using namespace Rcpp;


//============================ variance ========================================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec variance(const char * family,
                   const arma::vec & mu_vector) {
  arma::vec ans(mu_vector.n_elem);
  if(std::strcmp(family, "gaussian") == 0){
    ans.fill(1.0);
  }else if(std::strcmp(family, "binomial") == 0){
    ans = mu_vector % (1.0 - mu_vector);
  }else if(std::strcmp(family, "poisson") == 0){
    ans = mu_vector;
  }else if(std::strcmp(family, "Gamma") == 0){
    ans = pow(mu_vector, 2.0);
  }else if(std::strcmp(family, "inverse.gaussian") == 0){
    ans = pow(mu_vector, 3.0);
  }
  return(ans);
}
//==============================================================================


//============================ derivative variance wrt mean ====================
// [[Rcpp::export]]
arma::vec variancemu(const char * family,
                     const arma::vec & mu_vector) {
  arma::vec ans(mu_vector.n_elem);
  if(std::strcmp(family, "gaussian") == 0){
    ans.fill(0.0);
  }else if(std::strcmp(family, "poisson") == 0){
    ans.fill(1.0);
  }else if(std::strcmp(family, "binomial") == 0){
    ans = 1.0 - 2.0 * mu_vector;
  }else if(std::strcmp(family, "Gamma") == 0){
    ans = 2.0 * mu_vector;
  }else if(std::strcmp(family, "inverse.gaussian") == 0){
    ans = 3.0 * pow(mu_vector, 2.0);
  }
  return(ans);
}
//==============================================================================


//============================ second derivative variance wrt mean =============
// [[Rcpp::export]]
arma::vec variancemu2(const char * family,
                      const arma::vec & mu_vector) {
  arma::vec ans(mu_vector.n_elem);
  if(std::strcmp(family, "gaussian") == 0){
    ans.fill(0.0);
  }else if(std::strcmp(family, "poisson") == 0){
    ans.fill(0.0);
  }else if(std::strcmp(family, "binomial") == 0){
    ans.fill(-2.0);
  }else if(std::strcmp(family, "Gamma") == 0){
    ans.fill(2.0);
  }else if(std::strcmp(family, "inverse.gaussian") == 0){
    ans = 6.0 * mu_vector;
  }
  return(ans);
}
//==============================================================================
