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

#endif // UTILS_POLYNOMIAL_SELECTION