#include <stdlib.h>
#include <gmp.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"

// Operations on one polynomial

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

// Operations on two polynomials

void poly_prod(polynomial_mpz *res, polynomial_mpz f, polynomial_mpz g)
{
    mpz_t tmp;
    mpz_init_set_ui(tmp, 0);

    if ((f.degree == 0 && !mpz_cmp_ui(f.coeffs[0], 0)) || (g.degree == 0 && !mpz_cmp_ui(g.coeffs[0], 0)))
    {
        res->degree = 0;
        set_coeff(res, tmp, 0);
    }
    else
    {
        res->degree = f.degree + g.degree;

        for (size_t i = 0 ; i <= res->degree ; i++)
        {
            set_coeff(res, tmp, res->degree - i); // Reverse traversal so that we allocate memory once
        }

        for (size_t i = 0 ; i <= f.degree ; i++)
        {
            for (size_t j = 0 ; j <= g.degree ; j++)
            {
                mpz_addmul(res->coeffs[i+j], f.coeffs[i], g.coeffs[j]);
            }
        }
    }

    mpz_clear(tmp);
}

void poly_div(polynomial_mpz *res, polynomial_mpz f, polynomial_mpz g) // Returns remainder polynomial only, not quotient
{
    if (mpz_cmp_ui(g.coeffs[0], 1))
    {
        printf("Polynomial division where the divisor is not a monic polynomial.\n");
    }
    else if (f.degree >= g.degree)
    {
        if (g.degree) res->degree = g.degree - 1;
        else res->degree = 0;

        polynomial_mpz tmp_poly;
        init_poly_degree(&tmp_poly, f.degree);

        for (size_t i = 0 ; i <= f.degree ; i++)
        {
            set_coeff(&tmp_poly, f.coeffs[i], f.degree - i); // Same as poly prod, set in reversal order to allocate only once
        }

        mpz_t tmp, tmp2;
        mpz_inits(tmp, tmp2, NULL);

        mpz_set_ui(tmp2, 0);

        for (size_t i = 0 ; i <= f.degree - g.degree ; i++)
        {
            mpz_set(tmp, tmp_poly.coeffs[i]);

            set_coeff(&tmp_poly, tmp2, f.degree - i);

            for (size_t j = 0 ; j < g.degree ; j++)
            {
                mpz_submul(tmp_poly.coeffs[f.degree - i - j], tmp, g.coeffs[g.degree - j]);
            }
        }

        mpz_clears(tmp, tmp2, NULL);

        bool is_zero = true;
        for (size_t i = 0 ; i <= tmp_poly.degree ; i++)
        {
            if (mpz_cmp_ui(tmp_poly.coeffs[i], 0))
            {
                is_zero = false;
                break;
            }
        }

        if (is_zero)
        {
            reset_polynomial(res);
        }
        else
        {
            reduce_polynomial(&tmp_poly);

            for (size_t i = 0 ; i < g.degree ; i++)
            {
                set_coeff(res, tmp_poly.coeffs[i + f.degree - g.degree], g.degree - 1 - i);
            }
        }

        free_polynomial(&tmp_poly);

    } else
    {
        res->degree = f.degree;
        for (size_t i = 0 ; i <= f.degree ; i++)
        {
            set_coeff(res, f.coeffs[f.degree - i], f.degree - i); // Same as poly prod, set in reversal order to allocate only once
        }
    }
}