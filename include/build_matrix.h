#ifndef BUILD_MATRIX_H
#define BUILD_MATRIX_H

#include "dynamic_arrays.h"
#include "algebraic_base.h"
#include "quadratic_characters.h"
#include "NFS_relations.h"

void build_sparse_matrix(
    dyn_array_classic *sparse_matrix,
    nfs_relations *relations,
    dyn_array_classic *rational_primes,
    algebraic_base *algebraic_primes,
    quadratic_character_base *quad_char,
    unsigned long *divide_leading,
    unsigned long len_divide_leading
);

#endif // BUILD_MATRIX_H