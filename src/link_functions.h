#ifndef GEER_LINK_FUNCTIONS_H
#define GEER_LINK_FUNCTIONS_H

#include <RcppArmadillo.h>
#include "link_codes.h"
#include "family_codes.h"

namespace geer {

arma::vec inverse_link(LinkCode link_code,
                       const arma::vec& eta_vector);
arma::vec link_derivative_1(LinkCode link_code,
                            const arma::vec& eta_vector);
arma::vec link_derivative_2(LinkCode link_code,
                            const arma::vec& eta_vector);
arma::vec link_derivative_3(LinkCode link_code,
                            const arma::vec& eta_vector);
arma::vec link_derivative_1(const char* link,
                            const arma::vec& eta_vector);
// Validity checks on the linear predictor and the fitted means. They take the
// parsed codes so that the solvers can validate trial points cheaply.
bool is_valid_eta(LinkCode link_code,
                  const arma::vec& eta_vector);
bool is_valid_mu(FamilyCode family_code,
                 const arma::vec& mu_vector);

}  // namespace geer

#endif
