#define ARMA_WARN_LEVEL 1
#include "link_functions.h"
#include "family_utils.h"
#include "link_utils.h"
#include "utils.h"
#include <cfloat>
#include <cmath>


namespace {
inline bool all_finite(const Rcpp::NumericVector& x) {
  return Rcpp::is_true(Rcpp::all(Rcpp::is_finite(x)));
}
  inline bool all_positive_finite(const Rcpp::NumericVector& x) {
    return Rcpp::is_true(Rcpp::all(Rcpp::is_finite(x) & (x > 0.0)));
  }
  inline bool all_nonzero_finite(const Rcpp::NumericVector& x) {
    return Rcpp::is_true(Rcpp::all(Rcpp::is_finite(x) & (x != 0.0)));
  }
  inline bool all_open_unit_interval(const Rcpp::NumericVector& x) {
    return Rcpp::is_true(Rcpp::all(Rcpp::is_finite(x) & (x > 0.0) & (x < 1.0)));
  }
  inline arma::vec arma_logistic_mu(const arma::vec& eta) {
    const arma::vec eta_clipped = arma::clamp(eta, -30.0, 30.0);
    const arma::vec p = 1.0 / (1.0 + arma::exp(-eta_clipped));
    return arma::clamp(p, DBL_EPSILON, 1.0 - DBL_EPSILON);
  }
}


//============================ link inverse - arma (char*) ====================
arma::vec linkinv(const char* link,
                  const arma::vec& eta_vector) {
  return linkinv(parse_link(link), eta_vector);
}
//==============================================================================


//============================ link inverse - arma (enum) =====================
arma::vec linkinv(LinkCode lc, const arma::vec& eta) {
  const arma::uword n = eta.n_elem;
  switch (lc) {
  case LinkCode::logit:
    return arma_logistic_mu(eta);
  case LinkCode::probit: {
    const double thr = -R::qnorm(DBL_EPSILON, 0.0, 1.0, true, false);
    const arma::vec eta_clipped = arma::clamp(eta, -thr, thr);
    arma::vec ans(n);
    for (arma::uword i = 0; i < n; ++i)
      ans[i] = R::pnorm(eta_clipped[i], 0.0, 1.0, true, false);
    return ans;
  }
  case LinkCode::cauchit: {
    const double thr = -R::qcauchy(DBL_EPSILON, 0.0, 1.0, true, false);
    const arma::vec eta_clipped = arma::clamp(eta, -thr, thr);
    arma::vec ans(n);
    for (arma::uword i = 0; i < n; ++i)
      ans[i] = R::pcauchy(eta_clipped[i], 0.0, 1.0, true, false);
    return ans;
  }
  case LinkCode::cloglog: {
    const arma::vec eta_clipped = arma::clamp(eta, -arma::datum::inf, 700.0);
    return arma::clamp(1.0 - arma::exp(-arma::exp(eta_clipped)),
                       DBL_EPSILON, 1.0 - DBL_EPSILON);
  }
  case LinkCode::identity:
    return eta;
  case LinkCode::log:
    return arma::clamp(arma::exp(arma::clamp(eta, -arma::datum::inf, 700.0)),
                       DBL_EPSILON, arma::datum::inf);
  case LinkCode::sqrt:
    return eta % eta;
  case LinkCode::inverse_mu_squared:
    return 1.0 / arma::sqrt(eta);
  case LinkCode::inverse:
    return 1.0 / eta;
  }
  Rcpp::stop("Unsupported link.");
}
//==============================================================================



//============================ mu eta - first derivative - arma (char*) ========
arma::vec mueta(const char* link,
                const arma::vec& eta_vector) {
  return mueta(parse_link(link), eta_vector);
}
//==============================================================================


//============================ mu eta - first derivative - arma (enum) =========
arma::vec mueta(LinkCode lc, const arma::vec& eta) {
  const arma::uword n = eta.n_elem;
  switch (lc) {
  case LinkCode::logit: {
    const arma::vec mu = arma_logistic_mu(eta);
    return arma::clamp(mu % (1.0 - mu), DBL_EPSILON, arma::datum::inf);
  }
  case LinkCode::probit: {
    arma::vec ans(n);
    for (arma::uword i = 0; i < n; ++i)
      ans[i] = std::max(R::dnorm(eta[i], 0.0, 1.0, false), DBL_EPSILON);
    return ans;
  }
  case LinkCode::cauchit: {
    arma::vec ans(n);
    for (arma::uword i = 0; i < n; ++i)
      ans[i] = std::max(R::dcauchy(eta[i], 0.0, 1.0, false), DBL_EPSILON);
    return ans;
  }
  case LinkCode::cloglog: {
    const arma::vec eta_clipped = arma::clamp(eta, -arma::datum::inf, 700.0);
    return arma::clamp(arma::exp(eta_clipped - arma::exp(eta_clipped)),
                       DBL_EPSILON, arma::datum::inf);
  }
  case LinkCode::identity:
    return arma::ones<arma::vec>(n);
  case LinkCode::log:
    return arma::clamp(arma::exp(arma::clamp(eta, -arma::datum::inf, 700.0)),
                       DBL_EPSILON, arma::datum::inf);
  case LinkCode::sqrt:
    return arma::clamp(2.0 * eta, DBL_EPSILON, arma::datum::inf);
  case LinkCode::inverse_mu_squared:
    return -0.5 / arma::pow(eta, 1.5);
  case LinkCode::inverse:
    return -1.0 / arma::square(eta);
  }
  Rcpp::stop("Unsupported link.");
}
//==============================================================================



//============================ mu eta - second derivative - arma (char*) =======
arma::vec mueta2(const char* link,
                 const arma::vec& eta_vector) {
  return mueta2(parse_link(link), eta_vector);
}
//==============================================================================


//============================ mu eta - second derivative - arma (enum) ========
arma::vec mueta2(LinkCode lc, const arma::vec& eta) {
  const arma::uword n = eta.n_elem;
  switch (lc) {
  case LinkCode::logit: {
    const arma::vec mu = arma_logistic_mu(eta);
    const arma::vec me = arma::clamp(mu % (1.0 - mu), DBL_EPSILON, arma::datum::inf);
    return (1.0 - 2.0 * mu) % me;
  }
  case LinkCode::probit: {
    const arma::vec me = mueta(LinkCode::probit, eta);
    return -eta % me;
  }
  case LinkCode::cauchit: {
    const arma::vec me = mueta(LinkCode::cauchit, eta);
    return -2.0 * (eta / (arma::square(eta) + 1.0)) % me;
  }
  case LinkCode::cloglog: {
    const arma::vec eta_clipped = arma::clamp(eta, -arma::datum::inf, 700.0);
    const arma::vec me = mueta(LinkCode::cloglog, eta);
    return me % (1.0 - arma::exp(eta_clipped));
  }
  case LinkCode::identity:
    return arma::zeros<arma::vec>(n);
  case LinkCode::log:
    return arma::clamp(arma::exp(arma::clamp(eta, -arma::datum::inf, 700.0)),
                       DBL_EPSILON, arma::datum::inf);
  case LinkCode::sqrt:
    return arma::vec(n, arma::fill::value(2.0));
  case LinkCode::inverse_mu_squared:
    return 0.75 / arma::pow(eta, 2.5);
  case LinkCode::inverse:
    return 2.0 / arma::pow(eta, 3.0);
  }
  Rcpp::stop("Unsupported link.");
}
//==============================================================================



//============================ mu eta - third derivative - arma (char*) ========
arma::vec mueta3(const char* link,
                 const arma::vec& eta_vector) {
  return mueta3(parse_link(link), eta_vector);
}
//==============================================================================


//============================ mu eta - third derivative - arma (enum) =========
arma::vec mueta3(LinkCode lc, const arma::vec& eta) {
  const arma::uword n = eta.n_elem;
  switch (lc) {
  case LinkCode::logit: {
    const arma::vec mu = arma_logistic_mu(eta);
    const arma::vec me = arma::clamp(mu % (1.0 - mu), DBL_EPSILON, arma::datum::inf);
    return me % (1.0 - 6.0 * mu + 6.0 * arma::square(mu));
  }
  case LinkCode::probit: {
    const arma::vec me = mueta(LinkCode::probit, eta);
    return me % (arma::square(eta) - 1.0);
  }
  case LinkCode::cauchit: {
    const arma::vec me = mueta(LinkCode::cauchit, eta);
    return ((6.0 * arma::square(eta) - 2.0) /
            arma::square(arma::square(eta) + 1.0)) % me;
  }
  case LinkCode::cloglog: {
    const arma::vec eta_clipped = arma::clamp(eta, -arma::datum::inf, 350.0);
    const arma::vec exp_eta = arma::exp(eta_clipped);
    const arma::vec exp_2eta = arma::exp(2.0 * eta_clipped);
    const arma::vec me = mueta(LinkCode::cloglog, eta);
    return me % (1.0 - 3.0 * exp_eta + exp_2eta);
  }
  case LinkCode::identity:
  case LinkCode::sqrt:
    return arma::zeros<arma::vec>(n);
  case LinkCode::log:
    return arma::clamp(arma::exp(arma::clamp(eta, -arma::datum::inf, 700.0)),
                       DBL_EPSILON, arma::datum::inf);
  case LinkCode::inverse_mu_squared:
    return -1.875 / arma::pow(eta, 3.5);
  case LinkCode::inverse:
    return -6.0 / arma::pow(eta, 4.0);
  }
  Rcpp::stop("Unsupported link.");
}
//==============================================================================


//============================ valid eta =======================================
bool valideta(const char* link,
              const Rcpp::NumericVector& eta_vector) {
  switch (parse_link(link)) {
  case LinkCode::logit:
  case LinkCode::probit:
  case LinkCode::cauchit:
  case LinkCode::cloglog:
  case LinkCode::identity:
  case LinkCode::log:
    return all_finite(eta_vector);
  case LinkCode::sqrt:
  case LinkCode::inverse_mu_squared:
    return all_positive_finite(eta_vector);
  case LinkCode::inverse:
    return all_nonzero_finite(eta_vector);
  }
  Rcpp::stop("Unsupported link.");
}
//==============================================================================


//============================ valid mu ========================================
bool validmu(const char* family,
             const Rcpp::NumericVector& mu_vector) {
  switch (parse_family(family)) {
  case FamilyCode::gaussian:
    return all_finite(mu_vector);
  case FamilyCode::binomial:
    return all_open_unit_interval(mu_vector);
  case FamilyCode::poisson:
  case FamilyCode::gamma:
  case FamilyCode::inverse_gaussian:
    return all_positive_finite(mu_vector);
  }
  Rcpp::stop("Unsupported family.");
}
//==============================================================================
