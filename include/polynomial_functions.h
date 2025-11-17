#ifndef POLYNOMIAL_FUNCTIONS_H
#define POLYNOMIAL_FUNCTIONS_H

#include <gmp.h>

#include "polynomial_structures.h"

void poly_derivative(polynomial_mpz *res, polynomial_mpz f);
void evaluate_homogeneous(mpz_t res, polynomial_mpz f, mpz_t x, mpz_t y);

#endif // POLYNOMIAL_FUNCTIONS_H