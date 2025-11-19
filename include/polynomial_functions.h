#ifndef POLYNOMIAL_FUNCTIONS_H
#define POLYNOMIAL_FUNCTIONS_H

#include <gmp.h>

#include "polynomial_structures.h"

// Operations on one polynomial

void poly_derivative(polynomial_mpz *res, polynomial_mpz f);
unsigned long evaluate_mod_p(polynomial_mpz f, unsigned long x, unsigned long p);
void evaluate_homogeneous(mpz_t res, polynomial_mpz f, mpz_t x, mpz_t y);
void basic_find_roots(polynomial_mpz f, dyn_array_classic *roots, unsigned long p);
void power_poly_mod(polynomial_mpz *res, polynomial_mpz poly, polynomial_mpz f, unsigned long p, unsigned long exponent);
bool irreducible(polynomial_mpz f, unsigned long p)

// Operations on two polynomials

void poly_prod(polynomial_mpz *res, polynomial_mpz f, polynomial_mpz g);
void poly_div(polynomial_mpz *res, polynomial_mpz f, polynomial_mpz g);
void poly_div_mod(polynomial_mpz *res, polynomial_mpz f, polynomial_mpz g, unsigned long p);
void gcd_poly_mod(polynomial_mpz *res, polynomial_mpz *f, polynomial_mpz *g, unsigned long p);

#endif // POLYNOMIAL_FUNCTIONS_H