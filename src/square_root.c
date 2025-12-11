#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>

#include "logs.h"
#include "dynamic_arrays.h"
#include "NFS_relations.h"
#include "polynomial_functions.h"
#include "square_root.h"
#include "utils.h"

void extract_rational_square_root(
    mpz_t rational_square_root,
    const nfs_relations * restrict relations,
    const bool * restrict kernel_vector,
    const mpz_t n,
    const dyn_array_classic * restrict rational_primes
)
{
    unsigned long *prime_factors = calloc(rational_primes->len + 1, sizeof(unsigned long));

    unsigned long bit_count;

    mpz_t tmp_mpz, current_prime;
    mpz_inits(tmp_mpz, current_prime, NULL);

    for (size_t i = 0 ; i < relations->len ; i++)
    {
        if (kernel_vector[i])
        {
            if (mpz_cmp_ui(relations->rels[i].rational_norm, 0) < 0) prime_factors[0]++;

            for (size_t j = 0 ; j < rational_primes->len ; j++)
            {
                mpz_set_ui(current_prime, rational_primes->start[j]);
                bit_count = mpz_remove(tmp_mpz, relations->rels[i].rational_norm, current_prime);

                prime_factors[j+1] += bit_count;
            }
        }
    }

    if ((prime_factors[0]>>1)&1) mpz_set_si(rational_square_root, -1);
    else mpz_set_ui(rational_square_root, 1);

    for (size_t i = 0 ; i < rational_primes->len ; i++)
    {
        mpz_set_ui(current_prime, rational_primes->start[i]);

        mpz_powm_ui(tmp_mpz, current_prime, prime_factors[i+1]>>1, n);

        mpz_mul(rational_square_root, rational_square_root, tmp_mpz);

        mpz_mod(rational_square_root, rational_square_root, n);
    }

    mpz_clears(tmp_mpz, current_prime, NULL);
}

void extract_algebraic_square_root(
    mpz_t algebraic_square_root,
    const polynomial_mpz f_x,
    polynomial_mpz algebraic_square,
    const mpz_t m0,
    const mpz_t m1,
    const mpz_t leading_coeff,
    const mpz_t coeff_bound,
    const unsigned long inert_prime,
    gmp_randstate_t state,
    FILE *logfile
)
{
    polynomial_mpz algebraic_root;
    init_poly(&algebraic_root);

    square_root_poly_mod(&algebraic_root, algebraic_square, f_x, inert_prime, state);

    log_msg(logfile, "Initial root computed, lifting...");

    Newton_lift(&algebraic_root, &algebraic_square, f_x, coeff_bound, inert_prime);

    log_msg(logfile, "Root lifted.");

    mpz_t tmp_mpz;
    mpz_init_set(tmp_mpz, leading_coeff);
    mpz_mul(tmp_mpz, tmp_mpz, m0);

    evaluate_homogeneous(algebraic_square_root, algebraic_root, tmp_mpz, m1);

    mpz_clear(tmp_mpz);
}

void Newton_lift(
    polynomial_mpz *algebraic_root,
    polynomial_mpz *algebraic_square,
    const polynomial_mpz f,
    const mpz_t bound,
    const unsigned long p
)
{
    polynomial_mpz root;
    init_poly(&root);
    copy_polynomial(&root, algebraic_root);

    mpz_t modulo;
    mpz_init_set_ui(modulo, p);

    mpz_t tmp_mpz, tmp_mpz2;
    mpz_init_set(tmp_mpz, bound);
    mpz_mul_2exp(tmp_mpz, tmp_mpz, 1);

    mpz_init(tmp_mpz2);

    polynomial_mpz tmp_poly, tmp_poly2;
    init_poly(&tmp_poly);
    init_poly(&tmp_poly2);

    while (mpz_cmp(modulo, tmp_mpz) < 0)
    {
        mpz_mul(modulo, modulo, modulo);

        poly_prod(&tmp_poly, root, root);
        poly_div_mod_mpz(&tmp_poly2, tmp_poly, f, modulo);

        poly_prod(&tmp_poly, tmp_poly2, *algebraic_square);
        poly_div_mod_mpz(&tmp_poly2, tmp_poly, f, modulo);

        for (size_t i = 0 ; i <= tmp_poly2.degree ; i++)
        {
            mpz_neg(tmp_poly2.coeffs[i], tmp_poly2.coeffs[i]);
            mpz_mod(tmp_poly2.coeffs[i], tmp_poly2.coeffs[i], modulo);
        }

        mpz_add_ui(tmp_poly2.coeffs[tmp_poly2.degree], tmp_poly2.coeffs[tmp_poly2.degree], 3);
        mpz_mod(tmp_poly2.coeffs[tmp_poly2.degree], tmp_poly2.coeffs[tmp_poly2.degree], modulo);

        poly_prod(&tmp_poly, tmp_poly2, root);
        poly_div_mod_mpz(&tmp_poly2, tmp_poly, f, modulo);

        mpz_set_ui(tmp_mpz2, 2);
        mpz_invert(tmp_mpz2, tmp_mpz2, modulo);

        copy_polynomial(&root, &tmp_poly2);

        for (size_t i = 0 ; i <= root.degree ; i++)
        {
            mpz_mul(root.coeffs[i], root.coeffs[i], tmp_mpz2);
            mpz_mod(root.coeffs[i], root.coeffs[i], modulo);
        }

        poly_prod(&tmp_poly, root, *algebraic_square);
        poly_div_mod_mpz(&tmp_poly2, tmp_poly, f, modulo);
        

        mpz_div_2exp(tmp_mpz2, modulo, 1);

        for (size_t i = 0 ; i <= tmp_poly2.degree ; i++)
        {
            if (mpz_cmp(tmp_poly2.coeffs[i], tmp_mpz2) > 0)
            {
                mpz_sub(tmp_poly2.coeffs[i], tmp_poly2.coeffs[i], modulo);
            }
        }
    }

    copy_polynomial(algebraic_root, &tmp_poly2);

    mpz_clears(tmp_mpz, tmp_mpz2, modulo, NULL);

    free_polynomial(&tmp_poly);
    free_polynomial(&tmp_poly2);
    free_polynomial(&root);
}