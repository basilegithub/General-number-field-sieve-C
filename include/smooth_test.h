#ifndef SMOOTH_TEST_H
#define SMOOTH_TEST_H

#include <gmp.h>

#include "dynamic_arrays.h"
#include "single_linked_list.h"
#include "polynomial_structures.h"
#include "NFS_relations.h"

void pollard_rho(const mpz_t m, mpz_t p1, mpz_t p2, gmp_randstate_t state);

unsigned long build_product_tree(
    const dyn_array * restrict reported,
    dyn_array * restrict tree_array,
    const mpz_t prod_primes,
    const mpz_t prod_primes_p1,
    mpz_t tmp
);

void build_remainder_tree(
    const dyn_array * restrict reported,
    dyn_array * restrict tree_array,
    const mpz_t prod_primes,
    const mpz_t prod_primes_p1,
    unsigned long tmp_long
);

void batch_smooth(
    nfs_relations * restrict smooth_candidates,
    dyn_array *reported,
    dyn_array *tree_array,
    const mpz_t prod_primes,
    const mpz_t prod_primes_p1,
    const mpz_t limit,
    const mpz_t limit_2,
    const size_t index,
    gmp_randstate_t state
);

void naive_smooth(
    nfs_relations * restrict smooth_candidates,
    const dyn_array_classic * restrict primes,
    const mpz_t limit,
    const mpz_t limit_2,
    const size_t index,
    gmp_randstate_t state
);


#endif // SMOOTH_TEST_H