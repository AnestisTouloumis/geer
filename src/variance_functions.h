#ifndef VARIANCEFUNCTIONS_H
#define VARIANCEFUNCTIONS_H


arma::vec variance(const char * family,
                   const arma::vec & mu_vector);

arma::vec variancemu(const char * family,
                     const arma::vec & mu_vector);

arma::vec variancemu2(const char * family,
                      const arma::vec & mu_vector);


#endif
