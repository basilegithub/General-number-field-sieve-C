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

#endif // SQUARE_ROOT_H