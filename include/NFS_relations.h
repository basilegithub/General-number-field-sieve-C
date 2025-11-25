#ifndef NFS_RELATIONS_H
#define NFS_RELATIONS_H

#include <gmp.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"

typedef struct
{
    polynomial_mpz poly_g;
    polynomial_mpz poly_f;
    mpz_t rational_norm;
    dyn_array_classic rational_large_primes;
    mpz_t algebraic_norm;

} nfs_relation; // Structure that holds all the required information of a nfs relation

#endif // NFS_RELATIONS_H