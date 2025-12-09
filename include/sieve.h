#ifndef SIEVE_H
#define SIEVE_H

#include <gmp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "algebraic_base.h"
#include "quadratic_characters.h"
#include "NFS_relations.h"

void sieve(
    nfs_relations * restrict smooth_candidates,
    const polynomial_mpz sieve_poly,
    const algebraic_base * restrict alg_base,
    const dyn_array_classic * restrict logs,
    const mpz_t leading_coeff,
    const mpz_t m0,
    const mpz_t m1,
    const size_t len_divide_leading,
    const unsigned long b,
    const unsigned long offset,
    const unsigned long sieve_len,
    unsigned short * restrict sieve_array
);

#endif // SIEVE_H