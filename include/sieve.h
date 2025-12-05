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
    nfs_relations *smooth_candidates,
    polynomial_mpz sieve_poly,
    algebraic_base *alg_base,
    dyn_array_classic *logs,
    mpz_t leading_coeff,
    mpz_t m0,
    mpz_t m1,
    size_t len_divide_leading,
    unsigned long b,
    unsigned long offset,
    unsigned long sieve_len,
    unsigned short *sieve_array
);

#endif // SIEVE_H