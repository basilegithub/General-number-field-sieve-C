#ifndef INIT_FUNCTIONS_H
#define INIT_FUNCTIONS_H

#include <gmp.h>

#include "dynamic_arrays.h"

unsigned int compute_degree(mpz_t n, mpf_t ln2, mpf_t e);
void initialize_params(gmp_randstate_t state, mpf_t ln10, mpf_t ln2, mpf_t e);
void compute_smooth_bound(mpz_t n, mpz_t smooth_bound, mpf_t ln2, mpf_t e);
void compute_logs(dyn_array_classic *logs, dyn_array_classic primes);

void compute_free_relations(
    nfs_relations * restrict relations,
    const algebraic_base * restrict alg_base,
    const unsigned long * restrict divide_leading,
    mpz_t leading_coeff,
    mpz_t m1,
    const size_t len_divide_leading,
    const unsigned long degree
);

#endif // INIT_FUNCTIONS_H