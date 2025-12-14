#ifndef UTILS_POLYNOMIAL_SELECTION
#define UTILS_POLYNOMIAL_SELECTION

#include <gmp.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "polynomial_functions.h"

void get_Lnorm(
    mpz_t res,
    const polynomial_mpz * restrict F,
    const unsigned long skew,
    const unsigned long smoothness_bound,
    const mpf_t ln2,
    const mpf_t e
);

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

#endif // UTILS_POLYNOMIAL_SELECTION