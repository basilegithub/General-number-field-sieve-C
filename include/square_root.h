#ifndef SQUARE_ROOT_H
#define SQUARE_ROOT_H

#include <gmp.h>
#include <stdio.h>

#include "dynamic_arrays.h"
#include "NFS_relations.h"

void extract_rational_square_root(
    mpz_t rational_square_root,
    const nfs_relations * restrict relations,
    const bool * restrict kernel_vector,
    const mpz_t n,
    const dyn_array_classic * restrict rational_primes
);

void extract_algebraic_square_root(
    mpz_t algebraic_square_root,
    const polynomial_mpz * restrict f_x,
    polynomial_mpz * restrict algebraic_square,
    const mpz_t m0,
    const mpz_t m1,
    const mpz_t leading_coeff,
    const mpz_t coeff_bound,
    const unsigned long inert_prime,
    gmp_randstate_t state,
    FILE *logfile
);

void Newton_lift(
    polynomial_mpz *algebraic_root,
    polynomial_mpz *algebraic_square,
    const polynomial_mpz * restrict f,
    const mpz_t bound,
    const unsigned long p
);

#endif // SQUARE_ROOT_H