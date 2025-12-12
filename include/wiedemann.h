#ifndef WIEDEMANN_H
#define WIEDEMANN_H

#include <gmp.h>
#include <stdbool.h>

#include "dynamic_arrays.h"

void poly_anul(mpz_t D, mpz_t B, const unsigned long m);

void find_kernel_vectors(
    dyn_array * restrict kernel_vectors,
    dyn_array_classic * restrict A,
    mpz_t minimal_polynomial_estimate,
    size_t * restrict block,
    const size_t block_size,
    const unsigned long n,
    const unsigned long limit
);

void wiedemann(
    dyn_array * restrict kernel_vectors,
    dyn_array_classic * restrict A,
    mpz_t minimal_polynomial_estimate,
    const size_t block_size,
    const unsigned long n,
    const unsigned long limit
);

#endif // WIEDEMANN_H