#define ARMA_WARN_LEVEL 1
#include "family_utils.h"
#include "variance_functions.h"


//============================ variance (enum) ==================================
arma::vec variance(FamilyCode fc, const arma::vec& mu_vector) {
  arma::vec ans(mu_vector.n_elem);
  switch (fc) {
  case FamilyCode::gaussian:
    ans.fill(1.0);
    break;
  case FamilyCode::binomial:
    ans = mu_vector % (1.0 - mu_vector);
    break;
  case FamilyCode::poisson:
    ans = mu_vector;
    break;
  case FamilyCode::gamma:
    ans = arma::square(mu_vector);
    break;
  case FamilyCode::inverse_gaussian:
    ans = mu_vector % mu_vector % mu_vector;
    break;
  }
  return ans;
}
//==============================================================================


//============================ variance (char*) =================================
arma::vec variance(const char* family, const arma::vec& mu_vector) {
  return variance(parse_family(family), mu_vector);
}
//==============================================================================


//============================ derivative variance wrt mean (enum) =============
arma::vec variancemu(FamilyCode fc, const arma::vec& mu_vector) {
  arma::vec ans(mu_vector.n_elem);
  switch (fc) {
  case FamilyCode::gaussian:
    ans.fill(0.0);
    break;
  case FamilyCode::poisson:
    ans.fill(1.0);
    break;
  case FamilyCode::binomial:
    ans = 1.0 - 2.0 * mu_vector;
    break;
  case FamilyCode::gamma:
    ans = 2.0 * mu_vector;
    break;
  case FamilyCode::inverse_gaussian:
    ans = 3.0 * arma::square(mu_vector);
    break;
  }
  return ans;
}
//==============================================================================


//============================ derivative variance wrt mean (char*) ============
arma::vec variancemu(const char* family, const arma::vec& mu_vector) {
  return variancemu(parse_family(family), mu_vector);
}
//==============================================================================


//============================ second derivative variance wrt mean (enum) ======
arma::vec variancemu2(FamilyCode fc, const arma::vec& mu_vector) {
  arma::vec ans(mu_vector.n_elem);
  switch (fc) {
  case FamilyCode::gaussian:
  case FamilyCode::poisson:
    ans.fill(0.0);
    break;
  case FamilyCode::binomial:
    ans.fill(-2.0);
    break;
  case FamilyCode::gamma:
    ans.fill(2.0);
    break;
  case FamilyCode::inverse_gaussian:
    ans = 6.0 * mu_vector;
    break;
  }
  return ans;
}
//==============================================================================


//============================ second derivative variance wrt mean (char*) =====
arma::vec variancemu2(const char* family, const arma::vec& mu_vector) {
  return variancemu2(parse_family(family), mu_vector);
}
//==============================================================================
