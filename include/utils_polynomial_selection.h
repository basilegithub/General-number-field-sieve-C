#ifndef UTILS_POLYNOMIAL_SELECTION
#define UTILS_POLYNOMIAL_SELECTION

#include <gmp.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "polynomial_functions.h"

void evaluate_legendre(
    const unsigned long n,
    const double x,
    double * restrict P,
    double * restrict dP
);

void nodes_and_weights_gauss(
    const unsigned long n,
    double * restrict x,
    double * restrict w
);

void get_Lnorm(
    mpz_t res,
    const polynomial_mpz * restrict F,
    const unsigned long skew,
    const unsigned long smoothness_bound,
    const mpf_t ln2,
    const mpf_t e
);

double get_Escore(
    const polynomial_mpz * restrict f,
    const polynomial_mpz * restrict g,
    const double alpha,
    const unsigned long smoothness_bound,
    const unsigned long x_limit,
    const unsigned long y_limit,
    dickman_table * restrict table,
    const mpf_t ln2,
    const mpf_t e
);

void minimize_Lnorm(
    polynomial_mpz * restrict f,
    unsigned long skew_factor,
    const unsigned long smoothness_bound,
    mpz_t m0,
    mpz_t m1,
    mpf_t ln2,
    mpf_t e
);



#endif // UTILS_POLYNOMIAL_SELECTION