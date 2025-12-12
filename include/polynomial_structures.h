#ifndef POLYNOMIAL_STRUCTURES_H
#define POLYNOMIAL_STRUCTURES_H

#include <gmp.h>

typedef struct
{
    mpz_t* coeffs;
    unsigned long degree;
} polynomial_mpz; // Polynomial of arbitrary size coefficients

// Functions definition

void init_poly(polynomial_mpz * restrict polynomial);

void init_poly_degree(polynomial_mpz * restrict polynomial, const unsigned long degree);

void reduce_polynomial(polynomial_mpz * restrict polynomial);

void reduce_polynomial_last(polynomial_mpz * restrict polynomial);

void set_coeff(polynomial_mpz * restrict polynomial, const mpz_t number, const unsigned long index);

void copy_polynomial(polynomial_mpz *polynomial1, const polynomial_mpz *polynomial2);

void reset_polynomial(polynomial_mpz * restrict polynomial);

void free_polynomial(polynomial_mpz * restrict polynomial);

void print_polynomial(const polynomial_mpz * restrict polynomial);

#endif // POLYNOMIAL_STRUCTURES_H