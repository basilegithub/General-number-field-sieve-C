#ifndef SQUARE_ROOT_H
#define SQUARE_ROOT_H

#include <gmp.h>

#include "dynamic_arrays.h"
#include "NFS_relations.h"

void extract_rational_square_root(
    mpz_t rational_square_root,
    nfs_relations *relations,
    bool *kernel_vector,
    mpz_t n,
    dyn_array_classic *rational_primes
);

void extract_algebraic_square_root(
    mpz_t algebraic_square_root,
    polynomial_mpz f_x,
    polynomial_mpz algebraic_square,
    mpz_t m0,
    mpz_t m1,
    mpz_t leading_coeff,
    mpz_t coeff_bound,
    unsigned long inert_prime,
    gmp_randstate_t state
);

#endif // SQUARE_ROOT_H