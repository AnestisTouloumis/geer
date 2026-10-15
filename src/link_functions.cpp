#include "link_functions.h"
#include <algorithm>
#include <cfloat>
#include <cmath>


namespace {
arma::vec inverse_logit(const arma::vec& eta_vector) {
  const arma::vec eta_clipped = arma::clamp(eta_vector, -30.0, 30.0);
  const arma::vec p = 1.0 / (1.0 + arma::exp(-eta_clipped));
  return arma::clamp(p, DBL_EPSILON, 1.0 - DBL_EPSILON);
}
}  // namespace


//============================ link inverse - arma (enum) ======================
arma::vec inverse_link(LinkCode link_code, const arma::vec& eta_vector) {
  const arma::uword n = eta_vector.n_elem;
  switch (link_code) {
  case LinkCode::logit:
    return inverse_logit(eta_vector);
  case LinkCode::probit: {
    const double thr = -R::qnorm(DBL_EPSILON, 0.0, 1.0, true, false);
    const arma::vec eta_clipped = arma::clamp(eta_vector, -thr, thr);
    arma::vec result(n);
    for (arma::uword i = 0; i < n; ++i) {
      result[i] = R::pnorm(eta_clipped[i], 0.0, 1.0, true, false);
    }
    return result;
  }
  case LinkCode::cauchit: {
    const double thr = -R::qcauchy(DBL_EPSILON, 0.0, 1.0, true, false);
    const arma::vec eta_clipped = arma::clamp(eta_vector, -thr, thr);
    arma::vec result(n);
    for (arma::uword i = 0; i < n; ++i) {
      result[i] = R::pcauchy(eta_clipped[i], 0.0, 1.0, true, false);
    }
    return result;
  }
  case LinkCode::cloglog: {
    const arma::vec eta_clipped = arma::clamp(eta_vector, -arma::datum::inf, 700.0);
    return arma::clamp(1.0 - arma::exp(-arma::exp(eta_clipped)),
                       DBL_EPSILON, 1.0 - DBL_EPSILON);
  }
  case LinkCode::identity:
    return eta_vector;
  case LinkCode::log:
    return arma::clamp(arma::exp(arma::clamp(eta_vector, -arma::datum::inf, 700.0)),
                       DBL_EPSILON, arma::datum::inf);
  case LinkCode::sqrt:
    return eta_vector % eta_vector;
  case LinkCode::inverse_mu_squared:
    return 1.0 / arma::sqrt(eta_vector);
  case LinkCode::inverse:
    return 1.0 / eta_vector;
  }
  Rcpp::stop("unsupported link");
}
//==============================================================================


//============================ mu eta_vector - first derivative - arma (char*) ========
arma::vec link_derivative_1(const char* link,
                            const arma::vec& eta_vector) {
  return link_derivative_1(parse_link(link), eta_vector);
}
//==============================================================================


//============================ mu eta_vector - first derivative - arma (enum) =========
arma::vec link_derivative_1(LinkCode link_code, const arma::vec& eta_vector) {
  const arma::uword n = eta_vector.n_elem;
  switch (link_code) {
  case LinkCode::logit: {
    const arma::vec mu = inverse_logit(eta_vector);
    return arma::clamp(mu % (1.0 - mu), DBL_EPSILON, arma::datum::inf);
  }
  case LinkCode::probit: {
    arma::vec result(n);
    for (arma::uword i = 0; i < n; ++i) {
      result[i] = std::max(R::dnorm(eta_vector[i], 0.0, 1.0, false), DBL_EPSILON);
    }
    return result;
  }
  case LinkCode::cauchit: {
    arma::vec result(n);
    for (arma::uword i = 0; i < n; ++i) {
      result[i] = std::max(R::dcauchy(eta_vector[i], 0.0, 1.0, false), DBL_EPSILON);
    }
    return result;
  }
  case LinkCode::cloglog: {
    const arma::vec eta_clipped = arma::clamp(eta_vector, -arma::datum::inf, 700.0);
    return arma::clamp(arma::exp(eta_clipped - arma::exp(eta_clipped)),
                       DBL_EPSILON, arma::datum::inf);
  }
  case LinkCode::identity:
    return arma::ones<arma::vec>(n);
  case LinkCode::log:
    return arma::clamp(arma::exp(arma::clamp(eta_vector, -arma::datum::inf, 700.0)),
                       DBL_EPSILON, arma::datum::inf);
  case LinkCode::sqrt:
    return arma::clamp(2.0 * eta_vector, DBL_EPSILON, arma::datum::inf);
  case LinkCode::inverse_mu_squared:
    return -0.5 / arma::pow(eta_vector, 1.5);
  case LinkCode::inverse:
    return -1.0 / arma::square(eta_vector);
  }
  Rcpp::stop("unsupported link");
}
//==============================================================================


//============================ mu eta_vector - second derivative - arma (enum) ========
arma::vec link_derivative_2(LinkCode link_code, const arma::vec& eta_vector) {
  const arma::uword n = eta_vector.n_elem;
  switch (link_code) {
  case LinkCode::logit: {
    const arma::vec mu = inverse_logit(eta_vector);
    const arma::vec me = arma::clamp(mu % (1.0 - mu), DBL_EPSILON, arma::datum::inf);
    return (1.0 - 2.0 * mu) % me;
  }
  case LinkCode::probit: {
    const arma::vec me = link_derivative_1(LinkCode::probit, eta_vector);
    return -eta_vector % me;
  }
  case LinkCode::cauchit: {
    const arma::vec me = link_derivative_1(LinkCode::cauchit, eta_vector);
    return -2.0 * (eta_vector / (arma::square(eta_vector) + 1.0)) % me;
  }
  case LinkCode::cloglog: {
    const arma::vec eta_clipped = arma::clamp(eta_vector, -arma::datum::inf, 700.0);
    const arma::vec me = link_derivative_1(LinkCode::cloglog, eta_vector);
    return me % (1.0 - arma::exp(eta_clipped));
  }
  case LinkCode::identity:
    return arma::zeros<arma::vec>(n);
  case LinkCode::log:
    return arma::clamp(arma::exp(arma::clamp(eta_vector, -arma::datum::inf, 700.0)),
                       DBL_EPSILON, arma::datum::inf);
  case LinkCode::sqrt:
    return arma::vec(n, arma::fill::value(2.0));
  case LinkCode::inverse_mu_squared:
    return 0.75 / arma::pow(eta_vector, 2.5);
  case LinkCode::inverse:
    return 2.0 / arma::pow(eta_vector, 3.0);
  }
  Rcpp::stop("unsupported link");
}
//==============================================================================


//============================ mu eta_vector - third derivative - arma (enum) =========
arma::vec link_derivative_3(LinkCode link_code, const arma::vec& eta_vector) {
  const arma::uword n = eta_vector.n_elem;
  switch (link_code) {
  case LinkCode::logit: {
    const arma::vec mu = inverse_logit(eta_vector);
    const arma::vec me = arma::clamp(mu % (1.0 - mu), DBL_EPSILON, arma::datum::inf);
    return me % (1.0 - 6.0 * mu + 6.0 * arma::square(mu));
  }
  case LinkCode::probit: {
    const arma::vec me = link_derivative_1(LinkCode::probit, eta_vector);
    return me % (arma::square(eta_vector) - 1.0);
  }
  case LinkCode::cauchit: {
    const arma::vec me = link_derivative_1(LinkCode::cauchit, eta_vector);
    return ((6.0 * arma::square(eta_vector) - 2.0) /
            arma::square(arma::square(eta_vector) + 1.0)) % me;
  }
  case LinkCode::cloglog: {
    const arma::vec eta_clipped = arma::clamp(eta_vector, -arma::datum::inf, 350.0);
    const arma::vec exp_eta = arma::exp(eta_clipped);
    const arma::vec exp_2eta = arma::exp(2.0 * eta_clipped);
    const arma::vec me = link_derivative_1(LinkCode::cloglog, eta_vector);
    return me % (1.0 - 3.0 * exp_eta + exp_2eta);
  }
  case LinkCode::identity:
  case LinkCode::sqrt:
    return arma::zeros<arma::vec>(n);
  case LinkCode::log:
    return arma::clamp(arma::exp(arma::clamp(eta_vector, -arma::datum::inf, 700.0)),
                       DBL_EPSILON, arma::datum::inf);
  case LinkCode::inverse_mu_squared:
    return -1.875 / arma::pow(eta_vector, 3.5);
  case LinkCode::inverse:
    return -6.0 / arma::pow(eta_vector, 4.0);
  }
  Rcpp::stop("unsupported link");
}
//==============================================================================


//============================ valid eta_vector - arma (enum) =========================
bool is_valid_eta(LinkCode link_code,
                  const arma::vec& eta_vector) {
  const double* x = eta_vector.memptr();
  const arma::uword n = eta_vector.n_elem;
  switch (link_code) {
  case LinkCode::logit:
  case LinkCode::probit:
  case LinkCode::cauchit:
  case LinkCode::cloglog:
  case LinkCode::identity:
  case LinkCode::log:
    for (arma::uword i = 0; i < n; ++i) {
      if (!R_FINITE(x[i])) return false;
    }
    return true;
  case LinkCode::sqrt:
  case LinkCode::inverse_mu_squared:
    for (arma::uword i = 0; i < n; ++i) {
      if (!(R_FINITE(x[i]) && x[i] > 0.0)) return false;
    }
    return true;
  case LinkCode::inverse:
    for (arma::uword i = 0; i < n; ++i) {
      if (!(R_FINITE(x[i]) && x[i] != 0.0)) return false;
    }
    return true;
  }
  Rcpp::stop("unsupported link");
}
//==============================================================================


//============================ valid mu - arma (enum) ==========================
bool is_valid_mu(FamilyCode family_code,
                 const arma::vec& mu_vector) {
  const double* x = mu_vector.memptr();
  const arma::uword n = mu_vector.n_elem;
  switch (family_code) {
  case FamilyCode::gaussian:
    for (arma::uword i = 0; i < n; ++i) {
      if (!R_FINITE(x[i])) return false;
    }
    return true;
  case FamilyCode::binomial:
    for (arma::uword i = 0; i < n; ++i) {
      if (!(R_FINITE(x[i]) && x[i] > 0.0 && x[i] < 1.0)) return false;
    }
    return true;
  case FamilyCode::poisson:
  case FamilyCode::gamma:
  case FamilyCode::inverse_gaussian:
    for (arma::uword i = 0; i < n; ++i) {
      if (!(R_FINITE(x[i]) && x[i] > 0.0)) return false;
    }
    return true;
  }
  Rcpp::stop("unsupported family");
}
//==============================================================================

