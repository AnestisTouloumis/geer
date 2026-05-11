#ifndef LINK_FUNCTIONS_H
#define LINK_FUNCTIONS_H

#include <RcppArmadillo.h>
#include "link_utils.h"

arma::vec linkinv(LinkCode lc,
                  const arma::vec& eta_vector);
arma::vec mueta(LinkCode lc,
                const arma::vec& eta_vector);
arma::vec mueta2(LinkCode lc,
                 const arma::vec& eta_vector);
arma::vec mueta3(LinkCode lc,
                 const arma::vec& eta_vector);
arma::vec mueta(const char* link,
                const arma::vec& eta_vector);
bool valideta(const char* link,
              const Rcpp::NumericVector& eta_vector);
bool validmu(const char* family,
             const Rcpp::NumericVector& mu_vector);

#endif
