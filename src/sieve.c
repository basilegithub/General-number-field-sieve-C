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

void sieve(
    polynomial_mpz sieve_poly,
    dyn_array_classic rat_base,
    algebraic_base alg_base,
    dyn_array_classic logs,
    mpz_t leading_coeff,
    mpz_t m0,
    mpz_t m1,
    unsigned long b,
    unsigned long offset,
    unsigned long sieve_len,
    unsigned short *sieve_array
)
{
    algebraic_base_prime *sieve_prime;
    sieve_prime = alg_base.start;

    unsigned long p, r;

    unsigned long shift, log, rat_root, alg_root;

    size_t i = 0;

    mpz_t tmp, invmod_m1;
    mpz_inits(tmp, invmod_m1, NULL);

    while (sieve_prime != NULL)
    {
        p = sieve_prime->prime;

        shift = sieve_len%p;

        log = logs.start[i];

        mpz_set(tmp, m1);
        mpz_mod_ui(tmp, tmp, p);

        if (mpz_cmp_ui(tmp, 0)) // Rational sieve
        {
            mpz_set_ui(tmp, p);
            mpz_invert(invmod_m1, m1, tmp);
            mpz_mul(tmp, m0, invmod_m1);
            mpz_mod_ui(tmp, tmp, p);
            rat_root = (shift + (mpz_get_ui(tmp)*b)%p)%p;

            for (size_t j = rat_root ; j < 2*sieve_len ; j++)
            {
                sieve_array[j] += log;
            }
        }

        // Algebraic sieve

        for (size_t j = 0 ; j < sieve_prime->roots.len ; j++)
        {
            r = sieve_prime->roots.start[i];
            alg_root = (shift + (b*r)%p)%p;

            for (size_t k = alg_root ; k < 2*sieve_len ; k++)
            {
                sieve_array[k] += log;
            }
        }

        i++;

        sieve_prime = sieve_prime->next;
    }

    mpz_clears(tmp, invmod_m1, NULL);

    // Check for smooth candidates

    // For now use naive polynomial evaluation

    mpz_t rational_eval, algebraic_eval;
    mpz_inits(rational_eval, algebraic_eval, NULL);

    signed long a = -sieve_len;

    

    for (size_t i = 0 ; i < 2*sieve_len ; i++)
    {

    }

    mpz_clears(rational_eval, algebraic_eval, NULL);
}