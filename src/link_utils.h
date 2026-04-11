#ifndef LINK_UTILS_H
#define LINK_UTILS_H

#include <RcppArmadillo.h>
#include <string_view>

enum class LinkCode {
  logit,
  probit,
  cauchit,
  cloglog,
  identity,
  log,
  sqrt,
  inverse_mu_squared,
  inverse
};

inline LinkCode parse_link(const char* link) {
  const std::string_view lk(link);
  if (lk == "logit")    return LinkCode::logit;
  if (lk == "probit")   return LinkCode::probit;
  if (lk == "cauchit")  return LinkCode::cauchit;
  if (lk == "cloglog")  return LinkCode::cloglog;
  if (lk == "identity") return LinkCode::identity;
  if (lk == "log")      return LinkCode::log;
  if (lk == "sqrt")     return LinkCode::sqrt;
  if (lk == "1/mu^2")   return LinkCode::inverse_mu_squared;
  if (lk == "inverse")  return LinkCode::inverse;
  Rcpp::stop("Unsupported link: %s", link);
}

#endif
