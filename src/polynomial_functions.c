#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "polynomial_functions.h"
#include "utils.h"

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

void evaluate_poly(mpz_t res, polynomial_mpz f, signed long x)
{
    if (!x)
    {
        mpz_set(res, f.coeffs[f.degree]);
    }
    else
    {
        mpz_t tmp, tmp_res;
        mpz_inits(tmp, tmp_res, NULL);

        mpz_set_ui(tmp, x);
        mpz_pow_ui(tmp, tmp, f.degree);

        mpz_set_ui(tmp_res, 0);

        for (size_t i = 0 ; i < f.degree ; i++)
        {
            mpz_addmul(tmp_res, f.coeffs[i], tmp);
            mpz_divexact_ui(tmp, tmp, x);
        }

        mpz_add(tmp_res, tmp_res, f.coeffs[f.degree]);

        mpz_set(res, tmp_res);

        mpz_clears(tmp, tmp_res, NULL);
    }
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
    mpz_set_ui(res, 0);

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

void power_poly_mod(polynomial_mpz *res, polynomial_mpz poly, polynomial_mpz f, unsigned long p, unsigned long exponent)
{
    if (exponent == 1)
    {
        copy_polynomial(res, &poly);
        return;
    }
    else if (!exponent)
    {
        reset_polynomial(res);

        mpz_t tmp;
        mpz_init_set_ui(tmp, 1);

        set_coeff(res, tmp, 0);

        mpz_clear(tmp);
    }
    else
    {
        polynomial_mpz tmp_poly, tmp_poly2;
        init_poly(&tmp_poly);
        init_poly(&tmp_poly2);

        power_poly_mod(&tmp_poly, poly, f, p, exponent>>1);

        poly_prod(&tmp_poly2, tmp_poly, tmp_poly);

        poly_div_mod(&tmp_poly2, tmp_poly2, f, p);

        if (exponent&1)
        {
            poly_prod(&tmp_poly2, tmp_poly2, poly);
            poly_div_mod(&tmp_poly2, tmp_poly2, f, p);
        }

        copy_polynomial(res, &tmp_poly2);

        free_polynomial(&tmp_poly);
        free_polynomial(&tmp_poly2);
    }
}

void find_roots(polynomial_mpz f, dyn_array_classic *roots, unsigned long p, gmp_randstate_t state)
{
    polynomial_mpz reduced_f;
    init_poly_degree(&reduced_f, f.degree);

    copy_polynomial(&reduced_f, &f);

    for (size_t i = 0 ; i <= f.degree ; i++)
    {
        mpz_mod_ui(reduced_f.coeffs[i], reduced_f.coeffs[i], p);
    }

    reduce_polynomial(&reduced_f); // reduced_f is f (mod p)

    polynomial_mpz g, tmp_poly;
    init_poly(&g);
    init_poly_degree(&tmp_poly, 1);

    mpz_t tmp;
    mpz_init_set_ui(tmp, 1);

    set_coeff(&tmp_poly, tmp, 1);

    power_poly_mod(&g, tmp_poly, reduced_f, p, p); // Compute g = x^p (mod p)

    if (g.degree == 0)
    {
        set_coeff(&g, tmp, 1);
        mpz_neg(g.coeffs[0], g.coeffs[0]);
    }
    else
    {
        mpz_sub_ui(g.coeffs[g.degree-1], g.coeffs[g.degree-1], 1);
    } // Computes g - x = x^p - x (contains all linear factors)

    gcd_poly_mod(&tmp_poly, &reduced_f, &g, p);

    if (!mpz_cmp_ui(tmp_poly.coeffs[tmp_poly.degree], 0))
    {
        append_classic(roots, 0);
        reduce_polynomial(&tmp_poly);
    }

    second_step_roots(tmp_poly, roots, p, state);

    mpz_clear(tmp);

    free_polynomial(&tmp_poly);
    free_polynomial(&g);
    free_polynomial(&reduced_f);
}

void second_step_roots(polynomial_mpz f, dyn_array_classic *roots, unsigned long p, gmp_randstate_t state)
{
    if (f.degree == 0) return; // f = a -> 0 root

    if (f.degree == 1) // f = ax + b -> 1 root
    {
        mpz_t tmp, tmp2;
        mpz_inits(tmp, tmp2, NULL);

        mpz_set_ui(tmp2, p);

        mpz_invert(tmp, f.coeffs[0], tmp2);
        mpz_mul(tmp, tmp, f.coeffs[1]);
        mpz_neg(tmp, tmp);
        mpz_mod_ui(tmp, tmp, p);

        append_classic(roots, mpz_get_ui(tmp));

        mpz_clears(tmp, tmp2, NULL);

        return;
    }

    if (f.degree == 2) // f = ax^2 + bx + c -> either 0 or 1 or 2 roots
    {
        mpz_t tmp, tmp2, tmp3;
        mpz_inits(tmp, tmp2, tmp3, NULL);

        mpz_mul(tmp, f.coeffs[1], f.coeffs[1]);
        mpz_mul(tmp2, f.coeffs[0], f.coeffs[2]);
        mpz_submul_ui(tmp, tmp2, 4);

        if (!mpz_cmp_ui(tmp, 0))
        {
            mpz_set_ui(tmp2, p);
            mpz_mul_ui(tmp, f.coeffs[0], 2);
            mpz_invert(tmp, tmp, tmp2);
            mpz_mul(tmp, tmp, f.coeffs[1]);
            mpz_neg(tmp, tmp);
            mpz_mod_ui(tmp, tmp, p);

            append_classic(roots, mpz_get_ui(tmp));
        }
        else if (my_legendre(tmp, p) > -1)
        {
            sqrt_mod(tmp, p, state);

            mpz_set_ui(tmp3, p);
            mpz_mul_ui(tmp2, f.coeffs[0], 2); // tmp2 = 2a
            mpz_invert(tmp2, tmp2, tmp3); // tmp2 = (2a)^(-1) (mod p)
            mpz_mul(tmp, tmp, tmp2);
            mpz_mod_ui(tmp, tmp, p);

            mpz_mul(tmp3, tmp2, f.coeffs[1]);
            mpz_sub(tmp3, tmp, tmp3);
            mpz_mod_ui(tmp3, tmp3, p);
            append_classic(roots, mpz_get_ui(tmp3));

            mpz_submul_ui(tmp3, tmp, 2);
            mpz_mod_ui(tmp3, tmp3, p);
            append_classic(roots, mpz_get_ui(tmp3));
        }

        mpz_clears(tmp, tmp2, tmp3, NULL);

        return;
    }

    // If degree of f is higher than 2

    polynomial_mpz h;
    init_poly_degree(&h, 0);

    mpz_t a, tmp;
    mpz_inits(a, tmp, NULL);

    while (!h.degree || poly_equal(&h, &f))
    {
        mpz_set_ui(tmp, p);

        mpz_urandomm(a, state, tmp);

        polynomial_mpz tmp_poly;
        init_poly_degree(&tmp_poly, 1);

        set_coeff(&tmp_poly, a, 0);

        mpz_set_ui(tmp, 1);
        set_coeff(&tmp_poly, tmp, 1); // tmp = x + a

        power_poly_mod(&h, tmp_poly, f, p, (p-1)>>1);

        reduce_polynomial(&h);
        mpz_sub_ui(h.coeffs[h.degree], h.coeffs[h.degree], 1);
        mpz_mod_ui(h.coeffs[h.degree], h.coeffs[h.degree], p);

        gcd_poly_mod(&tmp_poly, &h, &f, p);
        copy_polynomial(&h, &tmp_poly);

        free_polynomial(&tmp_poly);
    }

    second_step_roots(h, roots, p, state);

    mpz_clears(a, tmp, NULL);
}

bool irreducible(polynomial_mpz f, unsigned long p)
{
    polynomial_mpz g;
    init_poly(&g);

    mpz_t tmp;
    mpz_init_set_ui(tmp, 1);

    set_coeff(&g, tmp, 1);

    mpz_set_ui(tmp, 0);
    set_coeff(&g, tmp, 0); // Just to be sure that g has all the correct coefficients

    polynomial_mpz tmp_poly, res_gcd;
    init_poly(&tmp_poly);
    init_poly(&res_gcd);

    for (size_t i = 0 ; i <= f.degree>>1 ; i++)
    {
        reset_polynomial(&tmp_poly);

        power_poly_mod(&tmp_poly, g, g, p, p);

        copy_polynomial(&g, &tmp_poly);
        if (tmp_poly.degree == 0)
        {
            mpz_set_ui(tmp, p - 1);
            set_coeff(&tmp_poly, tmp, 1);
        }
        else
        {
            mpz_sub_ui(tmp_poly.coeffs[tmp_poly.degree-1], tmp_poly.coeffs[tmp_poly.degree-1], 1);
            mpz_mod_ui(tmp_poly.coeffs[tmp_poly.degree-1], tmp_poly.coeffs[tmp_poly.degree-1], p);
        }

        reset_polynomial(&res_gcd);

        gcd_poly_mod(&res_gcd, &f, &tmp_poly, p);

        if (res_gcd.degree)
        {
            free_polynomial(&tmp_poly);
            free_polynomial(&res_gcd);

            free_polynomial(&g);

            mpz_clear(tmp);

            return false;
        }

    }

    free_polynomial(&tmp_poly);
    free_polynomial(&res_gcd);

    free_polynomial(&g);

    mpz_clear(tmp);

    return true;
}

// Operations on two polynomials

bool poly_equal(polynomial_mpz *f, polynomial_mpz *g)
{
    if (f->degree != g->degree) return false;

    for (size_t i = 0 ; i < f->degree ; i++)
    {
        if (mpz_cmp(f->coeffs[i], g->coeffs[i])) return false;
    }

    return true;
}

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

        for (size_t i = 0 ; i <= res->degree ; i++)
        {
            set_coeff(res, tmp, f.degree + g.degree - i); // Reverse traversal so that we allocate memory once
        }

        res->degree = f.degree + g.degree;

        for (size_t i = 0 ; i <= f.degree ; i++)
        {
            for (size_t j = 0 ; j <= g.degree ; j++)
            {
                mpz_addmul(res->coeffs[i + j], f.coeffs[i], g.coeffs[j]);
            }
        }
    }

    mpz_clear(tmp);
}

void poly_div(polynomial_mpz *res, polynomial_mpz f, polynomial_mpz g) // Returns remainder polynomial only, not quotient
{
    if (mpz_cmp_ui(g.coeffs[0], 1)) // If polynomial is non monic, it is unlikely that the division is possible
    {
        printf("Polynomial division where the divisor is not a monic polynomial.\n");
    }
    else if (f.degree >= g.degree)
    {
        polynomial_mpz remainder;
        init_poly_degree(&remainder, f.degree);

        for (size_t i = 0 ; i <= f.degree ; i++)
        {
            set_coeff(&remainder, f.coeffs[i], f.degree - i); // Same as poly prod, set in reversal order to allocate only once
        }

        mpz_t tmp, tmp2;
        mpz_inits(tmp, tmp2, NULL);

        mpz_set_ui(tmp2, 0);

        for (size_t i = 0 ; i <= f.degree - g.degree ; i++)
        {
            if (mpz_cmp_ui(remainder.coeffs[i], 0))
            {
                mpz_set(tmp, remainder.coeffs[i]);

                set_coeff(&remainder, tmp2, f.degree - i);

                for (size_t j = 0 ; j < g.degree ; j++)
                {
                    mpz_submul(remainder.coeffs[f.degree - i - j], tmp, g.coeffs[g.degree - j]);
                }
            }
        }

        mpz_clears(tmp, tmp2, NULL);

        // Check if remainder is zero

        bool is_zero = true;
        for (size_t i = 0 ; i <= remainder.degree ; i++)
        {
            if (mpz_cmp_ui(remainder.coeffs[i], 0))
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
            reduce_polynomial(&remainder);

            for (size_t i = 0 ; i <= remainder.degree ; i++)
            {
                set_coeff(res, remainder.coeffs[i], remainder.degree - i);
            }
        }

        free_polynomial(&remainder);

    }
    else // Just copy the polynomial, as no operation is needed
    {
        res->degree = f.degree;
        for (size_t i = 0 ; i <= f.degree ; i++)
        {
            set_coeff(res, f.coeffs[i], f.degree - i); // Same as poly prod, set in reversal order to allocate only once
        }
    }
}

void poly_div_mod(polynomial_mpz *res, polynomial_mpz f, polynomial_mpz g, unsigned long p) // Returns remainder polynomial only, not quotient
{
    if (f.degree >= g.degree)
    {
        // Copy&reduce f and reduce g mod p
        polynomial_mpz remainder, reduced_g;
        init_poly_degree(&remainder, f.degree);
        init_poly_degree(&reduced_g, g.degree);

        mpz_t tmp, tmp2, tmp3, prime;
        mpz_inits(tmp, tmp2, tmp3, prime, NULL);

        for (size_t i = 0 ; i <= f.degree ; i++)
        {
            mpz_mod_ui(tmp, f.coeffs[i], p);
            set_coeff(&remainder, tmp, f.degree - i); // Same as poly prod, set in reversal order to allocate only once
        }

        for (size_t i = 0 ; i <= g.degree ; i++)
        {
            mpz_mod_ui(tmp, g.coeffs[i], p);
            set_coeff(&reduced_g, tmp, g.degree - i);
        }

        // Delete leading 0 in reduced g

        reduce_polynomial(&reduced_g);

        // Perform long division

        mpz_set_ui(tmp2, 0);

        mpz_set_ui(prime, p);
        mpz_neg(tmp3, reduced_g.coeffs[0]);
        mpz_invert(tmp3, tmp3, prime);

        for (size_t i = 0 ; i <= f.degree - g.degree ; i++)
        {
            if (mpz_cmp_ui(remainder.coeffs[i], 0))
            {
                mpz_mul(tmp, remainder.coeffs[i], tmp3);
                mpz_mod(tmp, tmp, prime);

                set_coeff(&remainder, tmp2, f.degree - i);

                for (size_t j = 0 ; j < g.degree ; j++)
                {
                    mpz_addmul(remainder.coeffs[f.degree - i - j], tmp, g.coeffs[g.degree - j]);
                    mpz_mod(remainder.coeffs[f.degree - i - j], remainder.coeffs[f.degree - i - j], prime);
                }
            }
        }

        mpz_clears(tmp, tmp2, tmp3, prime, NULL);

        // Check if remainder is zero

        bool is_zero = true;
        for (size_t i = 0 ; i <= remainder.degree ; i++)
        {
            if (mpz_cmp_ui(remainder.coeffs[i], 0))
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
            reduce_polynomial(&remainder);

            for (size_t i = 0 ; i <= remainder.degree ; i++)
            {
                set_coeff(res, remainder.coeffs[i], remainder.degree - i);
            }
        }

        free_polynomial(&remainder);

    }
    else // Just reduce divided polynomial mod p
    {
        res->degree = f.degree;

        mpz_t tmp;
        mpz_init(tmp);

        for (size_t i = 0 ; i <= f.degree ; i++)
        {
            mpz_mod_ui(tmp, f.coeffs[i], p);
            set_coeff(res, tmp, f.degree - i); // Same as poly prod, set in reversal order to allocate only once
        }

        mpz_clear(tmp);
    }
}

void gcd_poly_mod(polynomial_mpz *res, polynomial_mpz *f, polynomial_mpz *g, unsigned long p)
{
    polynomial_mpz p1, p2, p3;
    init_poly_degree(&p1, f->degree);
    init_poly_degree(&p2, g->degree);
    init_poly(&p3);

    copy_polynomial(&p1, f);
    copy_polynomial(&p2, g);

    bool is_zero = true;
    for (size_t i = 0 ; i <= p2.degree ; i++)
    {
        if (mpz_cmp_ui(p2.coeffs[i], 0))
        {
            is_zero = false;
            break;
        }
    }

    while (!is_zero)
    {
        poly_div_mod(&p3, p1, p2, p);
        copy_polynomial(&p1, &p2);
        copy_polynomial(&p2, &p3);

        is_zero = true;
        for (size_t i = 0 ; i <= p2.degree ; i++)
        {
            if (mpz_cmp_ui(p2.coeffs[i], 0))
            {
                is_zero = false;
                break;
            }
        }
    }

    if (mpz_cmp_ui(p1.coeffs[0], 1))
    {
        mpz_t tmp, prime;
        mpz_init(tmp);
        mpz_init_set_ui(prime, p);
        
        mpz_invert(tmp, p1.coeffs[0], prime);

        for (size_t i = 0 ; i <= p1.degree ; i++)
        {
            mpz_mul(p1.coeffs[i], p1.coeffs[i], tmp);
            mpz_mod(p1.coeffs[i], p1.coeffs[i], prime);
        }

        mpz_clears(tmp, prime, NULL);
    }

    copy_polynomial(res, &p1);

    free_polynomial(&p1);
    free_polynomial(&p2);
    free_polynomial(&p3);
}