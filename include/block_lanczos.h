#ifndef BLOCK_LANCZOS_H
#define BLOCK_LANCZOS_H

#include "dynamic_arrays.h"

unsigned int switch_indices(const size_t d, const size_t mask);

void multiply_d(size_t * output, const size_t *dense_matrix, const size_t d, const size_t N);

void multiply_d_inplace(size_t * restrict dense_matrix, const size_t d, const size_t N);

void extract_columns(size_t * restrict W_inv, size_t * restrict d, size_t *T, size_t N);

void solve(mpz_t * restrict matrix, mpz_t * restrict kernel, const size_t nb_rows, const size_t matrix_len);

void block_lanczos(dyn_array * restrict output, dyn_array_classic * sparse_matrix, const size_t nb_relations, const size_t block_size, const unsigned long index, FILE *logfile);

#endif // BLOCK_LANCZOS_H