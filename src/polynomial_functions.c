#include <gmp.h>

#include "polynomial_structures.h"

void poly_derivative(polynomial_mpz *res, polynomial_mpz f)
{
    mpz_t tmp;
    mpz_init(tmp);

    for (size_t i = 0 ; i < f.degree ; i++)
    {
        mpz_set_ui(tmp, f.degree - i);
        mpz_mul(tmp, tmp, f.coeffs[i]);
        set_coeff(res, tmp, i);
    }

    res->degree = f.degree - 1;

    mpz_clear(tmp);
}

void evaluate_homogeneous(mpz_t res, polynomial_mpz f, mpz_t x, mpz_t y)
{
    mpz_set(res, 0);

    mpz_t tmp, tmp2, tmp3;
    mpz_inits(tmp, tmp2, tmp3, NULL);

    mpz_pow_ui(tmp, x, f.degree);
    mpz_set_ui(tmp2, 1);

    for (size_t i = 0 ; i < f.degree ; i++)
    {
        mpz_mul(tmp3, tmp, tmp2);
        mpz_addmul(res, tmp3, f.coeffs[i]);

        mpz_divexact(tmp, tmp, x);
        mpz_mul(tmp2, tmp2, y);
    }

    mpz_mul(tmp3, tmp, tmp2);
    mpz_addmul(res, tmp3, f.coeffs[f.degree]);

    mpz_clears(tmp, tmp2, tmp3, NULL);
}