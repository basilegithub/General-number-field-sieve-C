#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "dynamic_arrays.h"
#include "linear_algebra_utils.h"

void multiply_size_t(const dyn_array_classic * restrict A, const unsigned long n, const unsigned long index, const size_t * restrict b, size_t * restrict res)
{
    size_t * tmp = calloc(n, sizeof(size_t));
    unsigned long i = 0;
    size_t tmp2 = 0;

    for (unsigned long k = 0 ; k < A->len ; k++)
    {
        if (A->start[k] == index)
        {
            tmp[i] = tmp2;
            i++;
            tmp2 = 0;
        }
        else
        {
            tmp2 ^= b[A->start[k]];
        }
    }
    for (unsigned long j = 0 ; j < i ; j++) res[j] = tmp[j];
    for (unsigned long j = i ; j < n ; j++) res[j] = 0;

    free(tmp);
}

void multiply_sparse(const dyn_array_classic * restrict A, const unsigned long dim_out, const unsigned long index, const size_t * restrict b, size_t * restrict res)
{
    size_t * restrict tmp = calloc(dim_out, sizeof(size_t));
    unsigned long i = 0;
    size_t tmp2 = 0;

    for (unsigned long k = 0 ; k < A->len ; k++)
    {
        if (A->start[k] == index)
        {
            tmp[i] = tmp2;
            i++;
            tmp2 = 0;
        }
        else tmp2 ^= b[A->start[k]];
    }
    for (unsigned long j = 0 ; j < i ; j++) res[j] = tmp[j];
    for (unsigned long j = i ; j < dim_out ; j++) res[j] = 0;

    free(tmp);
}

size_t dot_prod(const unsigned long n, const bool * restrict lbd, const size_t * restrict x)
{
    size_t tmp = 0;
    for (unsigned long i = 0 ; i < n ; i++)
    {
        if (lbd[i]) tmp ^= x[i];
    }
    return tmp;
}

void add_vectors(size_t * restrict output, const size_t * restrict vec_a, const size_t *restrict vec_b, const size_t N)
{
    for (size_t i = 0 ; i < N ; i++)
    {
        output[i] = vec_a[i] ^ vec_b[i];
    }
}

void identity(size_t * restrict output, const size_t N)
{
    for (size_t i = 0 ; i < N ; i++)
    {
        output[i] = (1<<(N-i-1));
    }
}

void concatenate(size_t * restrict output, const size_t * restrict matrix_A, const size_t *restrict matrix_B, const size_t N, const size_t dim_out)
{
    for (size_t i = 0 ; i < dim_out ; i++)
    {
        output[i] = ((size_t)matrix_A[i]<<N) | matrix_B[i];
    }
}

void dense_multiply(size_t * restrict output, const size_t * restrict matrix_A, const size_t * restrict matrix_B, const size_t len_A, const size_t len_B)
{
    memset(output, 0, len_A*sizeof(size_t));

    size_t tmp;

    for (size_t i = 0 ; i < len_A ; i++)
    {
        tmp = 0;
        for (size_t j = 0 ; j < len_B ; j++)
        {
            tmp ^= matrix_B[j] * ((matrix_A[i] >> (len_B-j-1))&1);
        }
        output[i] = tmp;
    }
}

void sparse_multiply_transpose(const dyn_array_classic * restrict sparse_matrix, const size_t * restrict vector, size_t * restrict output, const unsigned long limit, const unsigned long dim)
{
    memset(output, 0, dim * sizeof(size_t));

    size_t row_index = 0;

    for (size_t i = 0 ; i < sparse_matrix->len ; i++)
    {
        if (sparse_matrix->start[i] == limit) row_index++;
        else
        {
            output[sparse_matrix->start[i]] ^= vector[row_index];
        }
    }
}

void dense_multiply_transpose(size_t * restrict output, size_t * restrict matrix, size_t * restrict vector, const size_t dim1, const size_t dim2) // computes transpose(matrix) * vector
{
    memset(output, 0, dim2 * sizeof(size_t));

    for (size_t i = 0 ; i < dim2 ; i++)
    {
        for (size_t j = 0 ; j < dim1 ; j++)
        {
            output[i] ^= vector[j] * ((matrix[j]>>(dim2-i-1))&1);
        }
    }
}

void transpose_dense(mpz_t * restrict output, size_t * restrict matrix, const size_t dim1, const size_t dim2) // computes transpose(matrix)
{

    for (size_t i = 0 ; i < dim2 ; i++)
    {
        mpz_set_ui(output[i], 0);

        for (size_t j = 0 ; j < dim1 ; j++)
        {
            if ((matrix[j]>>(dim2-i-1))&1)
            {
                mpz_setbit(output[i], dim1-j-1);
            }
        }
    }
}

void bin_poly_prod(mpz_t res, mpz_t poly_a, const mpz_t poly_b)
{
    mpz_t tmp_poly;
    mpz_init_set_ui(tmp_poly, 0);

    mpz_t tmp;
    mpz_init_set(tmp, poly_a);
    mpz_set_ui(tmp, mpz_sizeinbase(tmp, 2) - 1);

    for (unsigned long i = mpz_get_ui(tmp) ; i > 0 ; i--)
    {
        mpz_div_2exp(tmp, poly_a, i);
        if (mpz_odd_p(tmp)) mpz_xor(tmp_poly, tmp_poly, poly_b);
        mpz_mul_2exp(tmp_poly, tmp_poly, 1);
    }

    if (mpz_odd_p(poly_a)) mpz_xor(tmp_poly, tmp_poly, poly_b);
    mpz_set(res, tmp_poly);

    mpz_clears(tmp, tmp_poly, NULL);
}

void bin_div_poly(mpz_t quotient, mpz_t remainder, const mpz_t poly_a, const mpz_t poly_b)
{
    mpz_t A, B, tmp;
    mpz_inits(A, B, tmp, NULL);

    mpz_set(A, poly_a);
    mpz_set(B, poly_b);

    mpz_set_ui(quotient, 0);
    mpz_set(remainder, poly_a);

    size_t degA = mpz_sizeinbase(A, 2)-1;
    size_t degB = mpz_sizeinbase(B, 2)-1;

    while (mpz_cmp_ui(remainder, 0) && degA >= degB)
    {
        size_t shift = degA - degB;

        mpz_set_ui(tmp, 1);
        mpz_mul_2exp(tmp, tmp, shift);
        mpz_xor(quotient, quotient, tmp);

        mpz_mul_2exp(tmp, B, shift);
        mpz_xor(remainder, remainder, tmp);

        degA = mpz_sizeinbase(remainder, 2)-1;
    }

    mpz_clears(A, B, tmp, NULL);
}

void bin_gcd_poly(mpz_t res, const mpz_t poly_a, const mpz_t poly_b)
{
    mpz_t a, b, q, r;
    mpz_inits(a, b, q, r, NULL);

    mpz_set(a, poly_a);
    mpz_set(b, poly_b);

    while (mpz_cmp_ui(b, 0))
    {
        bin_div_poly(q, r, a, b);
        mpz_set(a, b);
        mpz_set(b, r);
    }

    mpz_set(res, a);

    mpz_clears(a, b, q, r, NULL);
}