#ifndef GEER_FAMILY_CODES_H
#define GEER_FAMILY_CODES_H

#include <RcppArmadillo.h>
#include "utils.h"
#include <string_view>

enum class FamilyCode {
  gaussian,
  binomial,
  poisson,
  gamma,
  inverse_gaussian
};

inline FamilyCode parse_family(const char* family) {
  if (family == nullptr) {
    signal_input_error("'family' must not be NULL");
  }
  const std::string_view fam(family);
  if (fam == "gaussian")         return FamilyCode::gaussian;
  if (fam == "binomial")         return FamilyCode::binomial;
  if (fam == "poisson")          return FamilyCode::poisson;
  if (fam == "Gamma")            return FamilyCode::gamma;
  if (fam == "inverse.gaussian") return FamilyCode::inverse_gaussian;
  signal_input_error("'family' must be one of: gaussian, binomial, "
      "poisson, Gamma, inverse.gaussian; received '%s'", family);
}

#endif
