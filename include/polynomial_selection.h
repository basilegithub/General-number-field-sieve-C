#ifndef POLYNOMIAL_SELECTION_H
#define POLYNOMIAL_SELECTION_H

#include <gmp.h>

typedef struct
{
    polynomial_mpz poly;
    mpz_t m0;
    mpz_t m1;
} polynomial_ranking_element;

typedef struct
{
    polynomial_ranking_element *ranking;
    size_t size;
    size_t len;
} polynomial_ranking;

// Function definition

void init_ranking(
    polynomial_ranking * restrict ranking,
    size_t size
);

void append_element(
    polynomial_ranking * restrict ranking,
    polynomial_ranking_element * restrict element
);

void insert_element(
    polynomial_ranking * restrict ranking,
    polynomial_ranking_element * restrict element,
    size_t index
);

void free_ranking(
    polynomial_ranking * restrict ranking
);

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