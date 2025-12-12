#ifndef BUILD_MATRIX_H
#define BUILD_MATRIX_H

#include "dynamic_arrays.h"
#include "algebraic_base.h"
#include "quadratic_characters.h"
#include "NFS_relations.h"

void build_sparse_matrix(
    dyn_array_classic * restrict sparse_matrix,
    const nfs_relations * restrict relations,
    const dyn_array_classic * restrict rational_primes,
    const algebraic_base *restrict algebraic_primes,
    const quadratic_character_base * restrict quad_char,
    const unsigned long * restrict divide_leading,
    const unsigned long len_divide_leading
);

#endif // BUILD_MATRIX_H