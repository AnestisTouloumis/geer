#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
double plackett_bernoulli_p11_cpp(const double p1,
                                  const double p2,
                                  const double theta) {
  if (!R_finite(p1) || p1 < 0.0 || p1 > 1.0) {
    stop("p1 must be in [0, 1].");
  }
  if (!R_finite(p2) || p2 < 0.0 || p2 > 1.0) {
    stop("p2 must be in [0, 1].");
  }
  if (NumericVector::is_na(theta) || theta < 0.0) {
    stop("theta must be non-negative.");
  }

  const double tol = 1e-12;
  const double lower = std::max(0.0, p1 + p2 - 1.0);
  const double upper = std::min(p1, p2);

  // Degenerate margins
  if (p1 <= tol || p2 <= tol) return 0.0;
  if (p1 >= 1.0 - tol) return p2;
  if (p2 >= 1.0 - tol) return p1;

  // Limiting cases
  if (theta == 0.0) return lower;
  if (!R_finite(theta)) return upper;

  // Independence
  if (std::fabs(theta - 1.0) < tol) return p1 * p2;

  const double a = 1.0 + (p1 + p2) * (theta - 1.0);
  double disc = a * a - 4.0 * theta * (theta - 1.0) * p1 * p2;

  if (disc < 0.0 && disc > -tol) disc = 0.0;
  if (disc < 0.0) {
    stop("Negative discriminant encountered.");
  }

  double p11 = (a - std::sqrt(disc)) / (2.0 * (theta - 1.0));

  // Trim tiny numerical errors
  if (p11 < lower && p11 > lower - 1e-10) p11 = lower;
  if (p11 > upper && p11 < upper + 1e-10) p11 = upper;

  return p11;
}


// [[Rcpp::export]]
double plackett_bernoulli_p11_cpp1(const double p1,
                                  const double p2,
                                  const double theta,
                                  const double margin_tol = 1e-8,
                                  const double theta_tol  = 1e-12) {
  if (!R_finite(p1) || p1 < 0.0 || p1 > 1.0) {
    stop("p1 must be in [0, 1].");
  }
  if (!R_finite(p2) || p2 < 0.0 || p2 > 1.0) {
    stop("p2 must be in [0, 1].");
  }
  if (NumericVector::is_na(theta) || theta < 0.0) {
    stop("theta must be non-negative.");
  }
  if (margin_tol < 0.0) {
    stop("margin_tol must be non-negative.");
  }
  if (theta_tol < 0.0) {
    stop("theta_tol must be non-negative.");
  }

  const double boundary_dist =
    std::min(std::min(p1, 1.0 - p1), std::min(p2, 1.0 - p2));

  // Numerical fallback near the Bernoulli boundary
  if (boundary_dist < margin_tol) {
    return p1 * p2;
  }

  // Exact independence
  if (std::fabs(theta - 1.0) < theta_tol) {
    return p1 * p2;
  }

  // Limiting dependence cases
  const double lower = std::max(0.0, p1 + p2 - 1.0);
  const double upper = std::min(p1, p2);

  if (theta == 0.0) {
    return lower;
  }
  if (!R_finite(theta)) {
    return upper;
  }

  const double a = 1.0 + (p1 + p2) * (theta - 1.0);
  double disc = a * a - 4.0 * theta * (theta - 1.0) * p1 * p2;

  // Protect against tiny negative roundoff
  if (disc < 0.0 && disc > -1e-14) {
    disc = 0.0;
  }
  if (disc < 0.0) {
    stop("Negative discriminant encountered.");
  }

  double p11 = (a - std::sqrt(disc)) / (2.0 * (theta - 1.0));

  // Trim tiny numerical spill outside admissible interval
  if (p11 < lower && p11 > lower - 1e-10) p11 = lower;
  if (p11 > upper && p11 < upper + 1e-10) p11 = upper;

  return p11;
}

inline double plackett_bernoulli_p11_smooth_fast(const double p1,
                                                 const double p2,
                                                 const double theta,
                                                 const double margin_tol = 1e-8,
                                                 const double theta_tol  = 1e-12) {
  const double indep = p1 * p2;

  const double boundary_dist =
    std::min(std::min(p1, 1.0 - p1), std::min(p2, 1.0 - p2));

  double exact;

  if (std::fabs(theta - 1.0) < theta_tol) {
    exact = indep;
  } else {
    const double tm1 = theta - 1.0;
    const double a   = 1.0 + (p1 + p2) * tm1;
    double disc      = a * a - 4.0 * theta * tm1 * indep;

    if (disc < 0.0 && disc > -1e-14) disc = 0.0;

    exact = (a - std::sqrt(disc)) / (2.0 * tm1);
  }

  if (margin_tol <= 0.0) {
    return exact;
  }

  double x = boundary_dist / margin_tol;
  if (x <= 0.0) x = 0.0;
  else if (x >= 1.0) x = 1.0;

  const double w = x * x * (3.0 - 2.0 * x);

  return indep + w * (exact - indep);
}

// [[Rcpp::export]]
double plackett_bernoulli_p11_smooth_fast_cpp(const double p1,
                                              const double p2,
                                              const double theta,
                                              const double margin_tol = 1e-8,
                                              const double theta_tol  = 1e-12) {
  return plackett_bernoulli_p11_smooth_fast(p1, p2, theta,
                                            margin_tol, theta_tol);
}
