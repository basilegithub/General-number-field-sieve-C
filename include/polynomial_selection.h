#ifndef POLYNOMIAL_SELECTION_H
#define POLYNOMIAL_SELECTION_H

#include <gmp.h>

void basic_polynomial_selection(polynomial_mpz * restrict polynomial, const mpz_t n, mpz_t m0, mpz_t m1, const unsigned long d);

void Kleinjung_poly_selection(
    polynomial_mpz * restrict res_poly,
    mpz_t m0,
    mpz_t m1,
    const mpz_t n,
    const dyn_array_classic * restrict primes,
    const unsigned long nb_roots,
    const unsigned long prime_bounds,
    const unsigned long multiplier,
    const mpz_t M,
    const unsigned long d,
    const unsigned long nb_poly_coarse_eval,
    const unsigned long nb_poly_precise_eval,
    gmp_randstate_t state,
    FILE *logfile
);

#endif // POLYNOMIAL_SELECTION_H