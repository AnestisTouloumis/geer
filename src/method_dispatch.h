#ifndef METHOD_DISPATCH_H
#define METHOD_DISPATCH_H

#include <RcppArmadillo.h>
#include <cstring>

enum class MethodCode {
  gee = 0,
  brgee_naive,
  brgee_robust,
  brgee_empirical,
  pgee_jeffreys
};

inline MethodCode method_code(const char* method) {
  if (method == nullptr) {
    Rcpp::stop("method is NULL");
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
  Rcpp::stop("Unknown method: %s", method);
  return MethodCode::gee; // never reached
}

#endif
