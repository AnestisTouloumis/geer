#ifndef FAMILY_UTILS_H
#define FAMILY_UTILS_H

#include <RcppArmadillo.h>
#include <string_view>

enum class FamilyCode {
  gaussian,
  binomial,
  poisson,
  gamma,
  inverse_gaussian
};

inline FamilyCode parse_family(const char* family) {
  const std::string_view fam(family);
  if (fam == "gaussian")         return FamilyCode::gaussian;
  if (fam == "binomial")         return FamilyCode::binomial;
  if (fam == "poisson")          return FamilyCode::poisson;
  if (fam == "Gamma")            return FamilyCode::gamma;
  if (fam == "inverse.gaussian") return FamilyCode::inverse_gaussian;
  Rcpp::stop("Unsupported family: %s", family);
}

#endif
