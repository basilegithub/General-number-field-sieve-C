#ifndef EXTRACT_SOLUTION_H
#define EXTRACT_SOLUTION_H

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "NFS_relations.h"

void extract_solution(
    nfs_relations *relations,
    bool *kernel_vector,
    dyn_array_classic *rational_primes,
    polynomial_mpz f_x,
    polynomial_mpz f_prime_sq,
    mpz_t n,
    mpz_t m0,
    mpz_t m1,
    mpz_t f_prime_eval,
    unsigned long inert_prime,
    unsigned long max_a_size
);

#endif // EXTRACT_SOLUTION_H