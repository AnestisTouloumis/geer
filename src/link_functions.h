#ifndef LINK_FUNCTIONS_H
#define LINK_FUNCTIONS_H

#include <RcppArmadillo.h>
#include "link_codes.h"
#include "family_codes.h"

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
// Validity checks on the linear predictor and the fitted means. They take the
// parsed codes so that the solvers can validate trial points cheaply.
bool valideta(LinkCode lc,
              const arma::vec& eta_vector);
bool validmu(FamilyCode fc,
             const arma::vec& mu_vector);

#endif
