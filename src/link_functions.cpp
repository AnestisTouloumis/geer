#include "link_functions.h"
#include "family_utils.h"
#include "link_utils.h"
#include "utils.h"

#include <cfloat>

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

  inline Rcpp::NumericVector clamp_probabilities(const Rcpp::NumericVector& x) {
    return Rcpp::pmax(Rcpp::pmin(x, 1.0 - DBL_EPSILON), DBL_EPSILON);
  }

  inline Rcpp::NumericVector logistic_mu(const Rcpp::NumericVector& eta_vector) {
    const Rcpp::NumericVector eta_clipped =
      Rcpp::pmin(Rcpp::pmax(eta_vector, -30.0), 30.0);
    const Rcpp::NumericVector ans = 1.0 / (1.0 + exp(-eta_clipped));
    return clamp_probabilities(ans);
  }

}  // namespace


//============================ link function - rcpp ============================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector linkfun_rcpp(const char* link,
                                 const Rcpp::NumericVector& mu_vector) {
  switch (parse_link(link)) {
  case LinkCode::logit:
    return Rcpp::qlogis(mu_vector, 0.0, 1.0, true, false);

  case LinkCode::probit:
    return Rcpp::qnorm(mu_vector, 0.0, 1.0, true, false);

  case LinkCode::cauchit:
    return Rcpp::qcauchy(mu_vector, 0.0, 1.0, true, false);

  case LinkCode::cloglog:
    return log(-log(1.0 - mu_vector));

  case LinkCode::identity:
    return mu_vector;

  case LinkCode::log:
    return log(mu_vector);

  case LinkCode::sqrt:
    return sqrt(mu_vector);

  case LinkCode::inverse_mu_squared:
    return 1.0 / (mu_vector * mu_vector);

  case LinkCode::inverse:
    return 1.0 / mu_vector;
  }

  Rcpp::stop("Unsupported link.");
}
//==============================================================================


//============================  link function - arma ===========================
// [[Rcpp::export]]
arma::vec linkfun(const char* link,
                  const arma::vec& mu_vector) {
  return vec2arma(linkfun_rcpp(link, arma2vec(mu_vector)));
}
//==============================================================================


//============================ link inverse - rcpp =============================
// [[Rcpp::export]]
Rcpp::NumericVector linkinv_rcpp(const char* link,
                                 const Rcpp::NumericVector& eta_vector) {
  switch (parse_link(link)) {
  case LinkCode::logit:
    return logistic_mu(eta_vector);

  case LinkCode::probit: {
    const double threshold = -R::qnorm(DBL_EPSILON, 0.0, 1.0, true, false);
    const Rcpp::NumericVector eta_clipped =
      Rcpp::pmin(Rcpp::pmax(eta_vector, -threshold), threshold);
    return Rcpp::pnorm(eta_clipped, 0.0, 1.0, true, false);
  }

  case LinkCode::cauchit: {
    const double threshold = -R::qcauchy(DBL_EPSILON, 0.0, 1.0, true, false);
    const Rcpp::NumericVector eta_clipped =
      Rcpp::pmin(Rcpp::pmax(eta_vector, -threshold), threshold);
    return Rcpp::pcauchy(eta_clipped, 0.0, 1.0, true, false);
  }

  case LinkCode::cloglog: {
    const Rcpp::NumericVector eta_clipped = Rcpp::pmin(eta_vector, 700.0);
    const Rcpp::NumericVector ans = 1.0 - exp(-exp(eta_clipped));
    return clamp_probabilities(ans);
  }

  case LinkCode::identity:
    return eta_vector;

  case LinkCode::log:
    return Rcpp::pmax(exp(Rcpp::pmin(eta_vector, 700.0)), DBL_EPSILON);

  case LinkCode::sqrt:
    return eta_vector * eta_vector;

  case LinkCode::inverse_mu_squared:
    return 1.0 / sqrt(eta_vector);

  case LinkCode::inverse:
    return 1.0 / eta_vector;
  }

  Rcpp::stop("Unsupported link.");
}
//==============================================================================


//============================ link inverse - arma =============================
// [[Rcpp::export]]
arma::vec linkinv(const char* link,
                  const arma::vec& eta_vector) {
  return vec2arma(linkinv_rcpp(link, arma2vec(eta_vector)));
}
//==============================================================================


//============================ mu eta - first derivative - rcpp ================
// [[Rcpp::export]]
Rcpp::NumericVector mueta_rcpp(const char* link,
                               const Rcpp::NumericVector& eta_vector) {
  switch (parse_link(link)) {
  case LinkCode::logit: {
    const Rcpp::NumericVector mu_vector = logistic_mu(eta_vector);
    return Rcpp::pmax(mu_vector * (1.0 - mu_vector), DBL_EPSILON);
  }

  case LinkCode::probit:
    return Rcpp::pmax(Rcpp::dnorm(eta_vector, 0.0, 1.0, false), DBL_EPSILON);

  case LinkCode::cauchit:
    return Rcpp::pmax(Rcpp::dcauchy(eta_vector, 0.0, 1.0, false), DBL_EPSILON);

  case LinkCode::cloglog: {
    const Rcpp::NumericVector eta_clipped = Rcpp::pmin(eta_vector, 700.0);
    return Rcpp::pmax(exp(eta_clipped - exp(eta_clipped)), DBL_EPSILON);
  }

  case LinkCode::identity: {
    Rcpp::NumericVector ans(eta_vector.size());
    ans.fill(1.0);
    return ans;
  }

  case LinkCode::log:
    return Rcpp::pmax(exp(Rcpp::pmin(eta_vector, 700.0)), DBL_EPSILON);

  case LinkCode::sqrt:
    return 2.0 * eta_vector;

  case LinkCode::inverse_mu_squared:
    return -0.5 / pow(eta_vector, 1.5);

  case LinkCode::inverse:
    return -1.0 / pow(eta_vector, 2.0);
  }

  Rcpp::stop("Unsupported link.");
}
//==============================================================================


//============================ mu eta - first derivative - arma ================
// [[Rcpp::export]]
arma::vec mueta(const char* link,
                const arma::vec& eta_vector) {
  return vec2arma(mueta_rcpp(link, arma2vec(eta_vector)));
}
//==============================================================================


//============================ mu eta - second derivative - rcpp ===============
// [[Rcpp::export]]
Rcpp::NumericVector mueta2_rcpp(const char* link,
                                const Rcpp::NumericVector& eta_vector) {
  switch (parse_link(link)) {
  case LinkCode::logit: {
    const Rcpp::NumericVector mu_vector = logistic_mu(eta_vector);
    const Rcpp::NumericVector mu_eta_vector =
      Rcpp::pmax(mu_vector * (1.0 - mu_vector), DBL_EPSILON);
    return (1.0 - 2.0 * mu_vector) * mu_eta_vector;
  }

  case LinkCode::probit:
    return -eta_vector * mueta_rcpp("probit", eta_vector);

  case LinkCode::cauchit:
    return -2.0 * (eta_vector / (pow(eta_vector, 2.0) + 1.0)) *
      mueta_rcpp("cauchit", eta_vector);

  case LinkCode::cloglog: {
    const Rcpp::NumericVector eta_clipped = Rcpp::pmin(eta_vector, 700.0);
    return mueta_rcpp("cloglog", eta_vector) * (1.0 - exp(eta_clipped));
  }

  case LinkCode::identity: {
    Rcpp::NumericVector ans(eta_vector.size());
    ans.fill(0.0);
    return ans;
  }

  case LinkCode::log:
    return Rcpp::pmax(exp(Rcpp::pmin(eta_vector, 700.0)), DBL_EPSILON);

  case LinkCode::sqrt: {
    Rcpp::NumericVector ans(eta_vector.size());
    ans.fill(2.0);
    return ans;
  }

  case LinkCode::inverse_mu_squared:
    return 0.75 / pow(eta_vector, 2.5);

  case LinkCode::inverse:
    return 2.0 / pow(eta_vector, 3.0);
  }

  Rcpp::stop("Unsupported link.");
}
//==============================================================================


//============================ mu eta - second derivative - arma ===============
// [[Rcpp::export]]
arma::vec mueta2(const char* link,
                 const arma::vec& eta_vector) {
  return vec2arma(mueta2_rcpp(link, arma2vec(eta_vector)));
}
//==============================================================================


//============================ mu eta - third derivative - rcpp ================
// [[Rcpp::export]]
Rcpp::NumericVector mueta3_rcpp(const char* link,
                                const Rcpp::NumericVector& eta_vector) {
  switch (parse_link(link)) {
  case LinkCode::logit: {
    const Rcpp::NumericVector mu_vector = logistic_mu(eta_vector);
    const Rcpp::NumericVector mu_eta_vector =
      Rcpp::pmax(mu_vector * (1.0 - mu_vector), DBL_EPSILON);
    return mu_eta_vector * (1.0 - 6.0 * mu_vector + 6.0 * pow(mu_vector, 2.0));
  }

  case LinkCode::probit:
    return mueta_rcpp("probit", eta_vector) * (pow(eta_vector, 2.0) - 1.0);

  case LinkCode::cauchit:
    return ((6.0 * pow(eta_vector, 2.0) - 2.0) /
            pow(pow(eta_vector, 2.0) + 1.0, 2.0)) *
              mueta_rcpp("cauchit", eta_vector);

  case LinkCode::cloglog: {
    const Rcpp::NumericVector eta_clipped = Rcpp::pmin(eta_vector, 350.0);
    const Rcpp::NumericVector exp_eta_vector = exp(eta_clipped);
    const Rcpp::NumericVector exp_2eta_vector = exp(2.0 * eta_clipped);
    return mueta_rcpp("cloglog", eta_vector) *
      (1.0 - 3.0 * exp_eta_vector + exp_2eta_vector);
  }

  case LinkCode::identity: {
    Rcpp::NumericVector ans(eta_vector.size());
    ans.fill(0.0);
    return ans;
  }

  case LinkCode::log:
    return Rcpp::pmax(exp(Rcpp::pmin(eta_vector, 700.0)), DBL_EPSILON);

  case LinkCode::sqrt: {
    Rcpp::NumericVector ans(eta_vector.size());
    ans.fill(0.0);
    return ans;
  }

  case LinkCode::inverse_mu_squared:
    return -1.875 / pow(eta_vector, 3.5);

  case LinkCode::inverse:
    return -6.0 / pow(eta_vector, 4.0);
  }

  Rcpp::stop("Unsupported link.");
}
//==============================================================================


//============================ mu eta - third derivative - arma ================
// [[Rcpp::export]]
arma::vec mueta3(const char* link,
                 const arma::vec& eta_vector) {
  return vec2arma(mueta3_rcpp(link, arma2vec(eta_vector)));
}
//==============================================================================


//============================ valid eta =======================================
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
