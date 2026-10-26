#ifndef GEER_METHOD_CODES_H
#define GEER_METHOD_CODES_H

#include <RcppArmadillo.h>
#include "utils.h"
#include <cstring>

enum class MethodCode {
  gee = 0,
  brgee_naive,
  brgee_robust,
  brgee_empirical,
  pgee_jeffreys
};

inline MethodCode parse_method(const char* method) {
  if (method == nullptr) {
    signal_input_error("'method' must not be NULL");
  }
  switch (method[0]) {
  case 'g':
    if (std::strcmp(method, "gee") == 0) return MethodCode::gee;
    break;
  case 'b':
    if (std::strcmp(method, "brgee-naive") == 0) return MethodCode::brgee_naive;
    if (std::strcmp(method, "brgee-robust") == 0) return MethodCode::brgee_robust;
    if (std::strcmp(method, "brgee-empirical") == 0) return MethodCode::brgee_empirical;
    break;
  case 'p':
    if (std::strcmp(method, "pgee-jeffreys") == 0) return MethodCode::pgee_jeffreys;
    break;
  }
  signal_input_error("'method' must be one of: gee, brgee-naive, "
      "brgee-robust, brgee-empirical, pgee-jeffreys; "
      "received '%s'", method);
}

#endif
