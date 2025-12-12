#ifndef POLYNOMIAL_SELECTION_H
#define POLYNOMIAL_SELECTION_H

#include <gmp.h>

void basic_polynomial_selection(polynomial_mpz * restrict polynomial, const mpz_t n, mpz_t m0, mpz_t m1, const unsigned long d);

#endif // POLYNOMIAL_SELECTION_H