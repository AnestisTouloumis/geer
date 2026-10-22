#ifndef GEER_NUISANCE_ESTIMATION_OR_H
#define GEER_NUISANCE_ESTIMATION_OR_H


#include <RcppArmadillo.h>

arma::vec get_marginalized_odds_ratios_or(const arma::vec& response_vector,
                                                 const arma::vec& id_vector,
                                                 const arma::vec& repeated_vector,
                                                 const arma::vec& weights_vector,
                                                 const double adding_constant,
                                                 const char* or_structure);

#endif
