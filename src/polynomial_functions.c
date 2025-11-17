#include <gmp.h>

#include "dynamic_arrays.h"
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

unsigned long evaluate_mod_p(polynomial_mpz f, unsigned long x, unsigned long p)
{
    mpz_t tmp, tmp_res;
    mpz_inits(tmp, tmp_res, NULL);

    if (!x)
    {
        mpz_set(tmp_res, f.coeffs[f.degree]);
        mpz_mod_ui(tmp_res, tmp_res, p);
    }
    else
    {
        mpz_set_ui(tmp, x);
        mpz_pow_ui(tmp, tmp, f.degree);

        mpz_set_ui(tmp_res, 0);

        for (size_t i = 0 ; i < f.degree ; i++)
        {
            mpz_addmul(tmp_res, f.coeffs[i], tmp);
            mpz_divexact_ui(tmp, tmp, x);
            mpz_mod_ui(tmp_res, tmp_res, p);
        }

        mpz_add(tmp_res, tmp_res, f.coeffs[f.degree]);
        mpz_mod_ui(tmp_res, tmp_res, p);
    }

    unsigned long res = mpz_get_ui(tmp_res);

    mpz_clears(tmp, tmp_res, NULL);

    return res;
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

void basic_find_roots(polynomial_mpz f, dyn_array_classic *roots, unsigned long p)
{
    roots->len = 0;

    for (unsigned long i = 0 ; i < p ; i++)
    {
        unsigned long eval = evaluate_mod_p(f, i, p);
        if (!eval)
        {
            append_classic(roots, i);
            if (roots->len == f.degree) return;
        }
    }
}