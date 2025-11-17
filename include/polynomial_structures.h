#ifndef POLYNOMIAL_STRUCTURES_H
#define POLYNOMIAL_STRUCTURES_H

#include <gmp.h>

typedef struct
{
    mpz_t* coeffs;
    unsigned long degree;
} polynomial_mpz; // Polynomial of arbitrary size coefficients

// Functions definition

void init_poly(polynomial_mpz *polynomial);
void init_poly_degree(polynomial_mpz *polynomial, unsigned long degree);
void set_coeff(polynomial_mpz *polynomial, mpz_t number, unsigned long index);
void copy_polynomial(polynomial_mpz *polynomial1, polynomial_mpz *polynomial2);
void poly_derivative(polynomial_mpz *res, polynomial_mpz *f);
void evaluate_homogeneous(mpz_t res, polynomial_mpz f, mpz_t x, mpz_t y)
void free_polynomial(polynomial_mpz *polynomial);
void print_polynomial(polynomial_mpz *polynomial);

#endif // POLYNOMIAL_STRUCTURES_H