#ifndef GEER_NUISANCE_ESTIMATION_OR_H
#define GEER_NUISANCE_ESTIMATION_OR_H

// Estimation of the marginalized odds ratios (the nuisance parameters of the
// odds-ratio GEE). The working covariance matrices built from them are
// declared in working_covariance_or.h.

#include <RcppArmadillo.h>

Rcpp::NumericVector get_marginalized_odds_ratios(const arma::vec& response_vector,
                                                 const arma::vec& id_vector,
                                                 const arma::vec& repeated_vector,
                                                 const arma::vec& weights_vector,
                                                 const double adding_constant,
                                                 const Rcpp::String& or_structure);

#endif
