#ifndef LINK_FUNCTIONS_H
#define LINK_FUNCTIONS_H

#include <RcppArmadillo.h>
#include "link_utils.h"

// ---- R-facing interface (char* -> parse -> compute) -------------------------
Rcpp::NumericVector linkfun_rcpp(const char* link,
                                 const Rcpp::NumericVector& mu_vector);
arma::vec linkfun(const char* link,
                  const arma::vec& mu_vector);

Rcpp::NumericVector linkinv_rcpp(const char* link,
                                 const Rcpp::NumericVector& eta_vector);
arma::vec linkinv(const char* link,
                  const arma::vec& eta_vector);

Rcpp::NumericVector mueta_rcpp(const char* link,
                               const Rcpp::NumericVector& eta_vector);
arma::vec mueta(const char* link,
                const arma::vec& eta_vector);

Rcpp::NumericVector mueta2_rcpp(const char* link,
                                const Rcpp::NumericVector& eta_vector);
arma::vec mueta2(const char* link,
                 const arma::vec& eta_vector);

Rcpp::NumericVector mueta3_rcpp(const char* link,
                                const Rcpp::NumericVector& eta_vector);
arma::vec mueta3(const char* link,
                 const arma::vec& eta_vector);

bool valideta(const char* link,
              const Rcpp::NumericVector& eta_vector);
bool validmu(const char* family,
             const Rcpp::NumericVector& mu_vector);

// ---- Hot-path enum overloads (no string parsing, no Rcpp round-trip) --------
arma::vec linkinv(LinkCode lc, const arma::vec& eta_vector);
arma::vec mueta(LinkCode lc,   const arma::vec& eta_vector);
arma::vec mueta2(LinkCode lc,  const arma::vec& eta_vector);
arma::vec mueta3(LinkCode lc,  const arma::vec& eta_vector);

#endif
