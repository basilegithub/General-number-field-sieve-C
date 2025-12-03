#include <gmp.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "algebraic_base.h"
#include "quadratic_characters.h"
#include "NFS_relations.h"
#include "logs.h"
#include "sieve.h"
#include "smooth_test.h"

void mono_cpu_sieve(
    nfs_relations *relations,
    polynomial_mpz f_x,
    polynomial_mpz g_x,
    dyn_array_classic rat_base,
    algebraic_base alg_base,
    size_t nb_Algebraic_pairs,
    quadratic_character_base quad_base,
    size_t nb_Quadratic_characters,
    mpz_t leading_coeff,
    mpz_t prod_primes,
    mpz_t m0,
    mpz_t m1,
    size_t sieve_len,
    mpz_t const1,
    mpz_t const2,
    unsigned long *divide_leading,
    mpz_t *pow_div,
    size_t len_divide_leading,
    dyn_array_classic logs,
    gmp_randstate_t state,
    FILE *logfile
)
{
    size_t required_relations = 3 + rat_base.len + nb_Algebraic_pairs + nb_Quadratic_characters + len_divide_leading;

    unsigned long offset = mpz_sizeinbase(const2, 2);

    log_msg(logfile, "Sieving started...");
    log_msg(logfile, "Need to collect at least %zu relations.", required_relations);
    log_blank_line(logfile);

    unsigned long b = 1;

    mpz_t gcd_b_cd;
    mpz_init(gcd_b_cd);

    mpz_t tmp;
    mpz_init(tmp);

    unsigned short *sieve_array = calloc(2*sieve_len, sizeof(unsigned short));

    while (relations->len < required_relations+20)
    {
        // Update contribution of gcd of b and c_d

        mpz_gcd_ui(gcd_b_cd, leading_coeff, b);
        unsigned long new_offset = offset + mpz_sizeinbase(gcd_b_cd, 2);

        // Sieve

        // Setup sieve

        nfs_relations smooth_candidates;
        init_relations(&smooth_candidates);

        polynomial_mpz sieve_poly;
        init_poly_degree(&sieve_poly, f_x.degree);

        copy_polynomial(&sieve_poly, &f_x);
        for (size_t i = 0 ; i <= sieve_poly.degree ; i++)
        {
            mpz_set_ui(tmp, b);
            mpz_pow_ui(tmp, tmp, i);
            mpz_mul(sieve_poly.coeffs[i], sieve_poly.coeffs[i], tmp);
        }

        memset(sieve_array, 0, 2*sieve_len*sizeof(unsigned short));

        // Perform sieve

        sieve(
            &smooth_candidates,
            sieve_poly,
            rat_base,
            alg_base,
            logs,
            leading_coeff,
            m0,
            m1,
            len_divide_leading,
            b,
            new_offset,
            sieve_len,
            sieve_array
        );

        // Verify smooth candidates

        naive_smooth(&smooth_candidates, rat_base, const1, const2, state);

        // Append true smooths to the collected relations

        for (size_t i = 0 ; i < smooth_candidates.len ; i++)
        {
            if (smooth_candidates.rels[i].nb_relations && smooth_candidates.rels[i].rational_large_primes.len == 0 && smooth_candidates.rels[i].algebraic_large_primes.start == NULL)
            {
                init_new_relation(relations, len_divide_leading);

                copy_polynomial(&relations->rels[relations->len-1].poly_g, &smooth_candidates.rels[i].poly_g);
                copy_polynomial(&relations->rels[relations->len-1].poly_f, &smooth_candidates.rels[i].poly_f);

                mpz_set(relations->rels[relations->len-1].algebraic_norm, smooth_candidates.rels[i].algebraic_norm);
                mpz_set(relations->rels[relations->len-1].rational_norm, smooth_candidates.rels[i].rational_norm);

                relations->rels[relations->len-1].nb_relations = 1;

                for (size_t j = 0 ; j < len_divide_leading ; j++)
                {
                    relations->rels[relations->len-1].divide_leading[j] = smooth_candidates.rels[i].divide_leading[j];
                }

                // relations->len++;
            }
        }

        // Handle partial relations

        // Increment b

        b++;

        // Cleanup

        clear_relations(&smooth_candidates);

        // Display evolution of relation collection process

        printf("\rb = %lu | %lu/(%lu+20) relations found",
            b,
            relations->len,
            required_relations
        );

        fflush(stdout);
    }

    printf("\n");

    mpz_clear(tmp);

    free(sieve_array);
}