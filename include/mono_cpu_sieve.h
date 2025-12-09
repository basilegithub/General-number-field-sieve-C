#ifndef MONO_CPU_SIEVE_H
#define MONO_CPU_SIEVE_H

#include <gmp.h>
#include <stdio.h>

#include "dynamic_arrays.h"

void mono_cpu_sieve(
    nfs_relations *relations,
    polynomial_mpz f_x,
    polynomial_mpz g_x,
    dyn_array_classic *rat_base,
    algebraic_base *alg_base,
    size_t nb_Algebraic_pairs,
    size_t nb_Quadratic_characters,
    mpz_t leading_coeff,
    mpz_t prod_primes,
    mpz_t m0,
    mpz_t m1,
    size_t sieve_len,
    mpz_t const1,
    mpz_t const2,
    unsigned long *divide_leading,
    mpz_t *pow_div,
    size_t len_divide_leading,
    dyn_array_classic *logs,
    gmp_randstate_t state,
    FILE *logfile,
    int flag_batch_smooth
);

#endif // MONO_CPU_SIEVE_H