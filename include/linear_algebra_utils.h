#ifndef LINEAR_ALGEBRA_UTILS_H
#define LINEAR_ALGEBRA_UTILS_H

#include <gmp.h>
#include <stdbool.h>

#include "dynamic_arrays.h"

void multiply_size_t(const dyn_array_classic * restrict A, const unsigned long n, const unsigned long index, const size_t * restrict b, size_t * restrict res);

void multiply_sparse(const dyn_array_classic * restrict A, const unsigned long dim_out, const unsigned long index, const size_t * restrict b, size_t * restrict res);

size_t dot_prod(const unsigned long n, const bool * restrict lbd, const size_t * restrict x);

void add_vectors(size_t * restrict output, const size_t * restrict vec_a, const size_t *restrict vec_b, const size_t N);

void identity(size_t * restrict output, const size_t N);

void concatenate(size_t * restrict output, const size_t * restrict matrix_A, const size_t *restrict matrix_B, const size_t N, const size_t dim_out);

void dense_multiply(size_t * restrict output, const size_t * restrict matrix_A, const size_t * restrict matrix_B, const size_t len_A, const size_t len_B);

void sparse_multiply_transpose(const dyn_array_classic * restrict sparse_matrix, const size_t * restrict vector, size_t * restrict output, const unsigned long limit, const unsigned long dim);

void dense_multiply_transpose(size_t * restrict output, size_t * restrict matrix, size_t * restrict vector, const size_t dim1, const size_t dim2);

void transpose_dense(mpz_t * restrict output, size_t * restrict matrix, const size_t dim1, const size_t dim2);

void bin_poly_prod(mpz_t res, mpz_t poly_a, const mpz_t poly_b);

void bin_div_poly(mpz_t quotient, mpz_t remainder, const mpz_t poly_a, const mpz_t poly_b);

void bin_gcd_poly(mpz_t res, const mpz_t poly_a, const mpz_t poly_b);

#endif // LINEAR_ALGEBRA_UTILS_H