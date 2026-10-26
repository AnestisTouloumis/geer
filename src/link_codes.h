#ifndef GEER_LINK_CODES_H
#define GEER_LINK_CODES_H

#include <RcppArmadillo.h>
#include "utils.h"
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
  if (link == nullptr) {
    signal_input_error("'link' must not be NULL");
  }
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
  signal_input_error("'link' must be one of: logit, probit, cauchit, "
      "cloglog, identity, log, sqrt, 1/mu^2, inverse; "
      "received '%s'", link);
}

#endif
