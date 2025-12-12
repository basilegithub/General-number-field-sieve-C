#ifndef POLYNOMIAL_FUNCTIONS_H
#define POLYNOMIAL_FUNCTIONS_H

#include <gmp.h>

#include "polynomial_structures.h"

// Operations on one polynomial

void poly_derivative(polynomial_mpz * restrict res, const polynomial_mpz * restrict f);

void evaluate_poly(mpz_t res, const polynomial_mpz * restrict f, const signed long x);

unsigned long evaluate_mod_p(const polynomial_mpz * restrict f, const unsigned long x, const unsigned long p);

void evaluate_homogeneous(mpz_t res, const polynomial_mpz * restrict f, const mpz_t x, const mpz_t y);

void basic_find_roots(const polynomial_mpz * restrict f, dyn_array_classic * restrict roots, const unsigned long p);

void power_poly_mod(polynomial_mpz * restrict res, const polynomial_mpz * restrict poly, const polynomial_mpz * restrict f, const unsigned long p, const unsigned long exponent);

void power_poly_mod_mpz(polynomial_mpz * restrict res, const polynomial_mpz * restrict poly, const polynomial_mpz * restrict f, const unsigned long p, const mpz_t exponent);

void find_roots(const polynomial_mpz * restrict f, dyn_array_classic * restrict roots, const unsigned long p, gmp_randstate_t state);

void second_step_roots(const polynomial_mpz * restrict f, dyn_array_classic * restrict roots, const unsigned long p, gmp_randstate_t state);

bool irreducible(const polynomial_mpz * restrict f, const unsigned long p);

int quadratic_residue(const polynomial_mpz * restrict poly, const polynomial_mpz * restrict f_x, const unsigned long p);

void square_root_poly_mod(polynomial_mpz * restrict res, const polynomial_mpz * restrict square, const polynomial_mpz * restrict f_x, const unsigned long p, gmp_randstate_t state);

// Operations on two polynomials

bool poly_equal(const polynomial_mpz * restrict f, const polynomial_mpz * restrict g);

void poly_prod(polynomial_mpz * restrict res, const polynomial_mpz *f, const polynomial_mpz *g);

void poly_div(polynomial_mpz * restrict res, const polynomial_mpz *f, const polynomial_mpz *g);

void poly_div_mod(polynomial_mpz * restrict res, const polynomial_mpz *f, const polynomial_mpz *g, const unsigned long p);

void poly_div_mod_mpz(polynomial_mpz * restrict res, const polynomial_mpz *f, const polynomial_mpz *g, const mpz_t p);

void quotient_poly_mod(polynomial_mpz * restrict res, const polynomial_mpz *f, const polynomial_mpz *g, const unsigned long p);

void gcd_poly_mod(polynomial_mpz * restrict res, const polynomial_mpz *f, const polynomial_mpz *g, const unsigned long p);

#endif // POLYNOMIAL_FUNCTIONS_H