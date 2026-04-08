#pragma once
#include <RcppArmadillo.h>

// Fast method dispatch (O(1) after first char check)
enum MethodCode {
  M_GEE = 0,
  M_BR_NAIVE,
  M_BR_ROBUST,
  M_BR_EMPIRICAL,
  M_PGEE_JEFFREYS
};
inline MethodCode method_code(const char* method) {
  if (method == nullptr) Rcpp::stop("method is NULL");
  switch (method[0]) {
  case 'g':
    if (std::strcmp(method, "gee") == 0) return M_GEE;
    break;
  case 'b':
    if (std::strcmp(method, "brgee-naive") == 0) return M_BR_NAIVE;
    if (std::strcmp(method, "brgee-robust") == 0) return M_BR_ROBUST;
    if (std::strcmp(method, "brgee-empirical") == 0) return M_BR_EMPIRICAL;
    break;
  case 'p':
    if (std::strcmp(method, "pgee-jeffreys") == 0) return M_PGEE_JEFFREYS;
    break;
  }
  Rcpp::stop("Unknown method: %s", method);
  return M_GEE; // never reached
}
