#ifndef LINEAR_ALGEBRA_UTILS_H
#define LINEAR_ALGEBRA_UTILS_H

#include <gmp.h>
#include <stdbool.h>

#include "dynamic_arrays.h"

void multiply_size_t(const dyn_array_classic A, const unsigned long n, const unsigned long index, size_t *b, size_t *res);
void multiply_sparse(const dyn_array_classic A, const unsigned long dim_out, const unsigned long index, const size_t *b, size_t *res);
size_t dot_prod(const unsigned long n, const bool *lbd, const size_t *x);
void add_vectors(size_t *output, const size_t * restrict vec_a, const size_t *restrict vec_b, const size_t N);
void identity(size_t *output, const size_t N);
void concatenate(size_t *output, const size_t * restrict matrix_A, const size_t *restrict matrix_B, const size_t N, const size_t dim_out);
void dense_multiply(size_t *output, const size_t *matrix_A, const size_t *matrix_B, const size_t len_A, const size_t len_B);
void sparse_multiply_transpose(const dyn_array_classic sparse_matrix, const size_t *vector, size_t *output, const unsigned long limit, const unsigned long dim);
void dense_multiply_transpose(size_t *output, size_t *matrix, size_t *vector, size_t dim1, size_t dim2);
void transpose_dense(mpz_t *output, size_t *matrix, size_t dim1, size_t dim2);
void bin_poly_prod(mpz_t res, mpz_t poly_a, const mpz_t poly_b);
void bin_div_poly(mpz_t quotient, mpz_t remainder, const mpz_t poly_a, const mpz_t poly_b);
void bin_gcd_poly(mpz_t res, mpz_t poly_a, mpz_t poly_b);

#endif // LINEAR_ALGEBRA_UTILS_H