#include <gmp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "algebraic_base.h"
#include "quadratic_characters.h"
#include "NFS_relations.h"
#include "polynomial_functions.h"
#include "utils.h"

void sieve(
    nfs_relations * restrict smooth_candidates,
    const polynomial_mpz sieve_poly,
    const algebraic_base * restrict alg_base,
    const dyn_array_classic * restrict logs,
    const mpz_t leading_coeff,
    const mpz_t m0,
    const mpz_t m1,
    const size_t len_divide_leading,
    const unsigned long b,
    const unsigned long offset,
    const unsigned long sieve_len,
    unsigned short * restrict sieve_array
)
{
    algebraic_base_prime *sieve_prime;
    sieve_prime = alg_base->start;

    unsigned long p, r;

    unsigned long shift, log, rat_root, alg_root;

    size_t i = 0;

    mpz_t tmp, invmod_m1;
    mpz_inits(tmp, invmod_m1, NULL);

    while (sieve_prime != NULL)
    {
        p = sieve_prime->prime;

        shift = sieve_len%p;

        log = logs->start[i];

        mpz_set(tmp, m1);
        mpz_mod_ui(tmp, tmp, p);

        if (mpz_cmp_ui(tmp, 0)) // Rational sieve
        {
            mpz_set_ui(tmp, p);
            mpz_invert(invmod_m1, m1, tmp);
            mpz_mul(tmp, m0, invmod_m1);
            mpz_mod_ui(tmp, tmp, p);
            rat_root = (shift + (mpz_get_ui(tmp)*b)%p)%p;

            #pragma omp simd
            for (unsigned short* ptr = sieve_array + rat_root ; ptr < sieve_array + 2*sieve_len ; ptr += p)
            {
                *ptr += log;
            }
        }

        // Algebraic sieve

        for (size_t j = 0 ; j < sieve_prime->roots.len ; j++)
        {
            r = sieve_prime->roots.start[j];
            alg_root = (shift + (b*r)%p)%p;

            #pragma omp simd
            for (unsigned short* ptr = sieve_array + alg_root ; ptr < sieve_array + 2*sieve_len ; ptr += p)
            {
                *ptr += log;
            }
        }

        i++;

        sieve_prime = sieve_prime->next;
    }

    mpz_clear(invmod_m1);

    // Check for smooth candidates

    // Use finite differences for polynomial evaluation

    mpz_t rational_eval, algebraic_eval, full_eval;
    mpz_inits(rational_eval, algebraic_eval, full_eval, NULL);

    signed long a = -sieve_len;

    mpz_mul_ui(rational_eval, m0, b);
    mpz_neg(rational_eval, rational_eval);
    mpz_submul_ui(rational_eval, m1, sieve_len);

    if (b&1)
    {
        mpz_t *algebraic_difs = calloc(sieve_poly.degree + 1, sizeof(mpz_t));
        for (size_t i = 0 ; i <= sieve_poly.degree ; i++)
        {
            mpz_init(algebraic_difs[i]);
            evaluate_poly(algebraic_difs[i], sieve_poly, a+i);
        }
        
        for (size_t q = 1 ; q <= sieve_poly.degree ; q++)
        {
            for (size_t k = sieve_poly.degree ; k > q-1 ; k--)
            {
                mpz_sub(algebraic_difs[k], algebraic_difs[k], algebraic_difs[k-1]);
            }
        }

        mpz_set(algebraic_eval, algebraic_difs[0]);

        for (size_t i = 0 ; i < 2*sieve_len ; i++)
        {
            if (a && !(gcd(abs(a), b) - 1) && mpz_cmp_ui(rational_eval, 0))
            {
                mpz_mul(full_eval, rational_eval, algebraic_eval);

                if (mpz_cmp_ui(full_eval, 0) && sieve_array[i] + offset > mpz_sizeinbase(full_eval, 2)-1)
                {
                    // Add smooth candidate
                    init_new_relation(smooth_candidates, len_divide_leading);

                    mpz_set_ui(tmp, b);
                    mpz_neg(tmp, tmp);
                    set_coeff(&smooth_candidates->rels[smooth_candidates->len-1].poly_g, tmp, 1); // p(x) = -b*x
                    set_coeff(&smooth_candidates->rels[smooth_candidates->len-1].poly_f, tmp, 1); // q(x) = -b*x

                    mpz_mul_si(tmp, leading_coeff, a);
                    set_coeff(&smooth_candidates->rels[smooth_candidates->len-1].poly_g, tmp, 0); // p(x) = c_d*a - b*x
                    mpz_set_si(tmp, a);
                    set_coeff(&smooth_candidates->rels[smooth_candidates->len-1].poly_f, tmp, 0); // q(x) = a - b*x

                    mpz_set(smooth_candidates->rels[smooth_candidates->len-1].rational_norm, rational_eval);

                    mpz_set(smooth_candidates->rels[smooth_candidates->len-1].algebraic_norm, algebraic_eval);

                    smooth_candidates->rels[smooth_candidates->len-1].nb_relations = 1;

                    for (size_t i = 0 ; i < len_divide_leading ; i++)
                    {
                        smooth_candidates->rels[smooth_candidates->len-1].divide_leading[i] = (bool)mpz_divisible_ui_p(leading_coeff, b);
                    }

                    // Large prime array and list are left with no large prime for now
                }
            }

            a++;

            mpz_add(algebraic_eval, algebraic_eval, algebraic_difs[1]);

            for (size_t q = 1 ; q < sieve_poly.degree ; q++)
            {
                mpz_add(algebraic_difs[q], algebraic_difs[q], algebraic_difs[q+1]);
            }

            mpz_add(rational_eval, rational_eval, m1);
        }

        for (size_t i = 0 ; i <= sieve_poly.degree ; i++) mpz_clear(algebraic_difs[i]);
        free(algebraic_difs);
    }
    else // If b is even, then a cannot be even -> skip them
    {
        size_t init = 0;

        if (!(a&1))
        {
            a++;
            init = 1;
            mpz_add(rational_eval, rational_eval, m1);
        }

        mpz_t *algebraic_difs = calloc(sieve_poly.degree + 1, sizeof(mpz_t));
        for (size_t i = 0 ; i <= sieve_poly.degree ; i++)
        {
            mpz_init(algebraic_difs[i]);
            evaluate_poly(algebraic_difs[i], sieve_poly, a + 2*i);
        }
        
        for (size_t q = 1 ; q <= sieve_poly.degree ; q++)
        {
            for (size_t k = sieve_poly.degree ; k > q-1 ; k--)
            {
                mpz_sub(algebraic_difs[k], algebraic_difs[k], algebraic_difs[k-1]);
            }
        }

        mpz_set(algebraic_eval, algebraic_difs[0]);

        for (size_t i = init ; i < 2*sieve_len ; i += 2)
        {
            if (a && !(gcd(abs(a), b) - 1) && mpz_cmp_ui(rational_eval, 0))
            {
                mpz_mul(full_eval, rational_eval, algebraic_eval);

                if (mpz_cmp_ui(full_eval, 0) && sieve_array[i] + offset > mpz_sizeinbase(full_eval, 2)-1)
                {
                    // Add smooth candidate
                    init_new_relation(smooth_candidates, len_divide_leading);

                    mpz_set_ui(tmp, b);
                    mpz_neg(tmp, tmp);
                    set_coeff(&smooth_candidates->rels[smooth_candidates->len-1].poly_g, tmp, 1); // p(x) = -b*x
                    set_coeff(&smooth_candidates->rels[smooth_candidates->len-1].poly_f, tmp, 1); // q(x) = -b*x

                    mpz_mul_si(tmp, leading_coeff, a);
                    set_coeff(&smooth_candidates->rels[smooth_candidates->len-1].poly_g, tmp, 0); // p(x) = c_d*a - b*x
                    mpz_set_si(tmp, a);
                    set_coeff(&smooth_candidates->rels[smooth_candidates->len-1].poly_f, tmp, 0); // q(x) = a - b*x

                    mpz_set(smooth_candidates->rels[smooth_candidates->len-1].rational_norm, rational_eval);

                    mpz_set(smooth_candidates->rels[smooth_candidates->len-1].algebraic_norm, algebraic_eval);

                    smooth_candidates->rels[smooth_candidates->len-1].nb_relations = 1;

                    for (size_t i = 0 ; i < len_divide_leading ; i++)
                    {
                        smooth_candidates->rels[smooth_candidates->len-1].divide_leading[i] = (bool)mpz_divisible_ui_p(leading_coeff, b);
                    }

                    // Large prime array and list are left with no large prime for now
                }
            }

            a += 2;

            mpz_add(algebraic_eval, algebraic_eval, algebraic_difs[1]);

            for (size_t q = 1 ; q < sieve_poly.degree ; q++)
            {
                mpz_add(algebraic_difs[q], algebraic_difs[q], algebraic_difs[q+1]);
            }

            mpz_addmul_ui(rational_eval, m1, 2);
        }
    }

    mpz_clears(rational_eval, algebraic_eval, full_eval, tmp, NULL);
}