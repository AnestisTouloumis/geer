#include "variance_functions.h"

namespace geer {


//============================ variance (enum) =================================
arma::vec variance_function(FamilyCode family_code, const arma::vec& mu_vector) {
  arma::vec result(mu_vector.n_elem);
  switch (family_code) {
  case FamilyCode::gaussian:
    result.fill(1.0);
    break;
  case FamilyCode::binomial:
    result = mu_vector % (1.0 - mu_vector);
    break;
  case FamilyCode::poisson:
    result = mu_vector;
    break;
  case FamilyCode::gamma:
    result = arma::square(mu_vector);
    break;
  case FamilyCode::inverse_gaussian:
    result = mu_vector % mu_vector % mu_vector;
    break;
  }
  return result;
}
//==============================================================================


//============================ variance (char*) ================================
arma::vec variance_function(const char* family, const arma::vec& mu_vector) {
  return variance_function(parse_family(family), mu_vector);
}
//==============================================================================


//============================ derivative variance wrt mean (enum) =============
arma::vec variance_derivative_1(FamilyCode family_code, const arma::vec& mu_vector) {
  arma::vec result(mu_vector.n_elem);
  switch (family_code) {
  case FamilyCode::gaussian:
    result.fill(0.0);
    break;
  case FamilyCode::poisson:
    result.fill(1.0);
    break;
  case FamilyCode::binomial:
    result = 1.0 - 2.0 * mu_vector;
    break;
  case FamilyCode::gamma:
    result = 2.0 * mu_vector;
    break;
  case FamilyCode::inverse_gaussian:
    result = 3.0 * arma::square(mu_vector);
    break;
  }
  return result;
}
//==============================================================================


//============================ second derivative variance wrt mean (enum) ======
arma::vec variance_derivative_2(FamilyCode family_code, const arma::vec& mu_vector) {
  arma::vec result(mu_vector.n_elem);
  switch (family_code) {
  case FamilyCode::gaussian:
  case FamilyCode::poisson:
    result.fill(0.0);
    break;
  case FamilyCode::binomial:
    result.fill(-2.0);
    break;
  case FamilyCode::gamma:
    result.fill(2.0);
    break;
  case FamilyCode::inverse_gaussian:
    result = 6.0 * mu_vector;
    break;
  }
  return result;
}
//==============================================================================

}  // namespace geer
