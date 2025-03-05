#include <RcppArmadillo.h>
using namespace Rcpp;


//============================ arma to vec =====================================
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::NumericVector arma2vec(const arma::vec & x) {
  Rcpp::NumericVector ans = Rcpp::NumericVector(x.begin(), x.end());
  return(ans);
}
//==============================================================================


//============================ vec to arma =====================================
// [[Rcpp::export]]
arma::vec vec2arma(const Rcpp::NumericVector & x) {
  arma::vec ans = Rcpp::as<arma::vec>(x);
  return(ans);
}
//==============================================================================


//============================ link function - rcpp ============================
// [[Rcpp::export]]
Rcpp::NumericVector linkfun_rcpp(const char * link,
                                 const Rcpp::NumericVector & mu_vector) {
  Rcpp::NumericVector ans;
  if(std::strcmp(link, "logit") == 0){
    ans = Rcpp::qlogis(mu_vector, 0.0, 1.0, true, false);
  }else if(std::strcmp(link, "probit") == 0){
    ans = Rcpp::qnorm(mu_vector, 0.0, 1.0, true, false);
  }else if(std::strcmp(link, "cauchit") == 0){
    ans = Rcpp::qcauchy(mu_vector, 0.0, 1.0, true, false);
  }else if(std::strcmp(link, "cloglog") == 0){
    ans = log(-log(1.0 - mu_vector));
  }else if(std::strcmp(link, "identity") == 0){
    ans = mu_vector;
  }else if(std::strcmp(link, "log") == 0){
    ans = log(mu_vector);
  }else if(std::strcmp(link, "sqrt") == 0){
    ans = sqrt(mu_vector);
  }else if(std::strcmp(link, "1/mu^2") == 0){
    ans = 1.0 / (mu_vector * mu_vector);
  }else if(std::strcmp(link, "inverse") == 0){
    ans = 1.0 / mu_vector;
  }else{
    Rcpp::stop("invalid link \n");
  }
  return(ans);
}
//==============================================================================


//============================  link function - arma ===========================
// [[Rcpp::export]]
arma::vec linkfun(const char * link,
                  const arma::vec & mu_vector) {
  arma::vec ans = vec2arma(linkfun_rcpp(link, arma2vec(mu_vector)));
  return(ans);
}
//==============================================================================


//============================ link inverse - rcpp =============================
// [[Rcpp::export]]
Rcpp::NumericVector linkinv_rcpp(const char * link,
                                 const Rcpp::NumericVector & eta_vector) {
  Rcpp::NumericVector ans;
  if(std::strcmp(link, "logit") == 0){
    Rcpp::NumericVector threshold(eta_vector.size());
    threshold.fill(30.0);
    Rcpp::NumericVector exp_eta_vector =
      Rcpp::ifelse(abs(eta_vector) < threshold, exp(eta_vector), eta_vector);
    exp_eta_vector =
      Rcpp::ifelse(eta_vector > threshold, 1/DBL_EPSILON, exp_eta_vector);
    exp_eta_vector =
      Rcpp::ifelse(eta_vector < -threshold, DBL_EPSILON, exp_eta_vector);
    ans = exp_eta_vector / (1.0 + exp_eta_vector);
  }else if(std::strcmp(link, "probit") == 0){
    Rcpp::NumericVector machine_eps(eta_vector.size());
    machine_eps.fill(DBL_EPSILON);
    Rcpp::NumericVector threshold = -Rcpp::qnorm(machine_eps, 0.0, 1.0, true, false);
    Rcpp::NumericVector eta_vector_new =
      Rcpp::pmin(Rcpp::pmax(eta_vector, -threshold), threshold);
    ans = Rcpp::pnorm(eta_vector_new, 0.0, 1.0, true, false);
  }else if(std::strcmp(link, "cauchit") == 0){
    Rcpp::NumericVector machine_eps(eta_vector.size());
    machine_eps.fill(DBL_EPSILON);
    Rcpp::NumericVector threshold = -Rcpp::qcauchy(machine_eps, 0.0, 1.0, true, false);
    Rcpp::NumericVector eta_vector_new =
      Rcpp::pmin(Rcpp::pmax(eta_vector, -threshold), threshold);
    ans = Rcpp::pcauchy(eta_vector_new, 0.0, 1.0, true, false);
  }else if(std::strcmp(link, "cloglog") == 0){
    Rcpp::NumericVector ans_old = 1.0 - exp(-exp(eta_vector));
    ans = Rcpp::pmax(Rcpp::pmin(ans_old, 1.0 - DBL_EPSILON), DBL_EPSILON);
  }else if(std::strcmp(link, "identity") == 0){
    ans = eta_vector;
  }else if(std::strcmp(link, "log") == 0){
    ans = Rcpp::pmax(exp(eta_vector), DBL_EPSILON);
  }else if(std::strcmp(link, "sqrt") == 0){
    ans = pow(eta_vector, 2.0);
  }else if(std::strcmp(link, "1/mu^2") == 0){
    ans = 1.0 / sqrt(eta_vector);
  }else if(std::strcmp(link, "inverse") == 0){
    ans = 1.0 / eta_vector;
  }else{
    Rcpp::stop("invalid link \n");
  }
  return(ans);
}
//==============================================================================


//============================ link inverse - arma =============================
// [[Rcpp::export]]
arma::vec linkinv(const char * link,
                  const Rcpp::NumericVector & eta_vector) {
  arma::vec ans = vec2arma(linkinv_rcpp(link, arma2vec(eta_vector)));
  return(ans);
}
//==============================================================================


//============================ mu eta - first derivative - rcpp ================
// [[Rcpp::export]]
Rcpp::NumericVector mueta_rcpp(const char * link,
                               const Rcpp::NumericVector & eta_vector) {
  Rcpp::NumericVector ans(eta_vector.size());
  if(std::strcmp(link, "logit") == 0){
    Rcpp::NumericVector threshold(eta_vector.size());
    threshold.fill(30.0);
    NumericVector one_plus_exp_vector = 1.0 + exp(eta_vector);
    ans = Rcpp::ifelse(Rcpp::abs(eta_vector) >= threshold,
                       DBL_EPSILON,
                       exp(eta_vector)/pow(one_plus_exp_vector, 2.0));
  }else if(std::strcmp(link, "probit") == 0){
    ans = Rcpp::pmax(Rcpp::dnorm(eta_vector, 0.0, 1.0, false), DBL_EPSILON);
  }else if(std::strcmp(link, "cauchit") == 0){
    ans = Rcpp::pmax(Rcpp::dcauchy(eta_vector, 0.0, 1.0, false), DBL_EPSILON);
  }else if(std::strcmp(link, "cloglog") == 0){
    Rcpp::NumericVector eta_vector_new = Rcpp::pmin(eta_vector, 700.0);
    ans = Rcpp::pmax(exp(eta_vector_new - exp(eta_vector_new)), DBL_EPSILON);
  }else if(std::strcmp(link, "identity") == 0){
    ans.fill(1.0);
  }else if(std::strcmp(link, "log") == 0){
    ans =  Rcpp::pmax(exp(eta_vector), DBL_EPSILON);
  }else if(std::strcmp(link, "sqrt") == 0){
    ans = 2.0 * eta_vector;
  }else if(std::strcmp(link, "1/mu^2") == 0){
    ans = -0.5 / pow(eta_vector, 1.5);
  }else if(std::strcmp(link, "inverse") == 0){
    ans = -1.0 / pow(eta_vector, 2.0);
  }else{
    Rcpp::stop("invalid link \n");
  }
  return(ans);
}
//==============================================================================


//============================ mu eta - first derivative - arma ================
// [[Rcpp::export]]
arma::vec mueta(const char * link,
                const arma::vec & eta_vector) {
  arma::vec ans = vec2arma(mueta_rcpp(link, arma2vec(eta_vector)));
  return(ans);
}
//==============================================================================


//============================ mu eta - second derivative - rcpp ===============
// [[Rcpp::export]]
Rcpp::NumericVector mueta2_rcpp(const char * link,
                                const Rcpp::NumericVector & eta_vector) {
  Rcpp::NumericVector ans(eta_vector.size());
  if(std::strcmp(link, "logit") == 0){
    ans = (1.0 - 2.0 * linkinv_rcpp("logit", eta_vector)) * mueta_rcpp("logit", eta_vector);
  }else if(std::strcmp(link, "probit") == 0){
    ans = -eta_vector * mueta_rcpp("probit", eta_vector);
  }else if(std::strcmp(link, "cauchit") == 0){
    ans = - 2.0 * (eta_vector/(pow(eta_vector, 2.0) + 1.0)) * mueta_rcpp("cauchit", eta_vector);
  }else if(std::strcmp(link, "cloglog") == 0){
    ans = mueta_rcpp("cloglog", eta_vector) * (1.0 - exp(eta_vector));
  }else if(std::strcmp(link, "identity") == 0){
    ans.fill(0);
  }else if(std::strcmp(link, "log") == 0){
    ans = Rcpp::pmax(exp(eta_vector), DBL_EPSILON);
  }else if(std::strcmp(link, "sqrt") == 0){
    ans.fill(2.0);
  }else if(std::strcmp(link, "1/mu^2") == 0){
    ans = 0.75 / pow(eta_vector, 2.5);
  }else if(std::strcmp(link, "inverse") == 0){
    ans = 2.0 / pow(eta_vector, 3.0);
  }else{
    Rcpp::stop("invalid link \n");
  }
  return(ans);
}
//==============================================================================


//============================ mu eta - second derivative - arma ===============
// [[Rcpp::export]]
arma::vec mueta2(const char * link,
                 const arma::vec & eta_vector){
  arma::vec ans = vec2arma(mueta2_rcpp(link, arma2vec(eta_vector)));
  return(ans);
}
//==============================================================================


//============================ mu eta - third derivative - rcpp ================
// [[Rcpp::export]]
Rcpp::NumericVector mueta3_rcpp(const char * link,
                                const Rcpp::NumericVector & eta_vector) {
  Rcpp::NumericVector ans(eta_vector.size());
  if(std::strcmp(link, "logit") == 0){
    Rcpp::NumericVector mueta_logis = mueta_rcpp("logit", eta_vector);
    ans = mueta_logis *
      ((exp(2.0) - 4.0) * mueta_logis +
      pow(1 - linkinv_rcpp("logit", eta_vector), 2.0));
  }else if(std::strcmp(link, "probit") == 0){
    ans = mueta_rcpp("probit", eta_vector) * (pow(eta_vector, 2.0) - 1.0);
  }else if(std::strcmp(link, "cauchit") == 0){
    ans = (6.0 * pow(eta_vector, 2.0) - 2.0) / pow(pow(eta_vector, 2.0) + 1.0, 2.0) *
      mueta_rcpp("cauchit", eta_vector);
  }else if(std::strcmp(link, "cloglog") == 0){
    Rcpp::NumericVector eta_vector_new = Rcpp::pmin(eta_vector, 700.0);
    Rcpp::NumericVector exp_eta_vector = exp(eta_vector_new);
    ans = mueta_rcpp("cloglog", eta_vector) *
      (exp(2) * exp_eta_vector + 1.0 - 3.0 * exp_eta_vector);
  }else if(std::strcmp(link, "identity") == 0){
    ans.fill(0);
  }else if(std::strcmp(link, "log") == 0){
    ans = Rcpp::pmax(exp(eta_vector), DBL_EPSILON);
  }else if(std::strcmp(link, "sqrt") == 0){
    ans.fill(0);
  }else if(std::strcmp(link, "1/mu^2") == 0){
    ans = - 1.875 / pow(eta_vector, 3.5);
  }else if(std::strcmp(link, "inverse") == 0){
    ans = - 6.0 / pow(eta_vector, 4.0);
  }else{
    Rcpp::stop("invalid link \n");
  }
  return(ans);
}
//==============================================================================


//============================ mu eta - third derivative - arma ================
// [[Rcpp::export]]
arma::vec mueta3(const char * link,
                 const arma::vec & eta_vector) {
  arma::vec ans = vec2arma(mueta3_rcpp(link, arma2vec(eta_vector)));
  return(ans);
}
//==============================================================================


//============================ valid eta =======================================
// [[Rcpp::export]]
Rcpp::LogicalVector valideta(const char * link,
                             const Rcpp::NumericVector & eta_vector) {
  Rcpp::LogicalVector ans;
  if(std::strcmp(link, "logit") == 0){
    ans = true;
  }else if(std::strcmp(link, "probit") == 0){
    ans = true;
  }else if(std::strcmp(link, "cauchit") == 0){
    ans = true;
  }else if(std::strcmp(link, "cloglog") == 0){
    ans = true;
  }else if(std::strcmp(link, "identity") == 0){
    ans = true;
  }else if(std::strcmp(link, "log") == 0){
    ans = true;
  }else if(std::strcmp(link, "sqrt") == 0){
    ans = is_true(all(Rcpp::is_finite(eta_vector) & (eta_vector > 0)));
  }else if(std::strcmp(link, "1/mu^2") == 0){
    ans = is_true(all(Rcpp::is_finite(eta_vector) & (eta_vector > 0)));
  }else if(std::strcmp(link, "inverse") == 0){
    ans = is_true(all(Rcpp::is_finite(eta_vector) & (eta_vector != 0)));
  }else{
    Rcpp::stop("invalid link \n");
  }
  return (ans);
}
//==============================================================================


//============================ valid mu ========================================
// [[Rcpp::export]]
Rcpp::LogicalVector validmu(const char * family,
                            const Rcpp::NumericVector & mu_vector) {
  Rcpp::LogicalVector ans;
  if(std::strcmp(family, "gaussian") == 0){
    ans = true;
  }else if(std::strcmp(family, "binomial") == 0){
    ans =
      is_true(all(Rcpp::is_finite(mu_vector) & (mu_vector > 0) & (mu_vector < 1)));
  }else if(std::strcmp(family, "poisson") == 0){
    ans = is_true(all(Rcpp::is_finite(mu_vector) & (mu_vector > 0)));
  }else if(std::strcmp(family, "Gamma") == 0){
    ans = is_true(all(Rcpp::is_finite(mu_vector) & (mu_vector > 0)));
  }else if(std::strcmp(family, "inverse.gaussian") == 0){
    ans = true;
  }else{
    Rcpp::stop("invalid family \n");
  }
  return(ans);
}
//==============================================================================
