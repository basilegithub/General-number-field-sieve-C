#ifndef EXTRACT_SOLUTION_H
#define EXTRACT_SOLUTION_H

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "NFS_relations.h"

void extract_solution(
    mpz_t factor,
    const nfs_relations * restrict relations,
    const bool * restrict kernel_vector,
    const dyn_array_classic * restrict rational_primes,
    const polynomial_mpz *f_x,
    const polynomial_mpz *f_prime_sq,
    const mpz_t leading_coeff,
    const mpz_t n,
    const mpz_t m0,
    const mpz_t m1,
    const mpz_t f_prime_eval,
    const unsigned long inert_prime,
    const unsigned long max_a_size,
    gmp_randstate_t state,
    FILE *logfile
);

#endif // EXTRACT_SOLUTION_H