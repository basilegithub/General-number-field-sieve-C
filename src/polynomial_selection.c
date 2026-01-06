#include <gmp.h>
#include <stdio.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "utils.h"
#include "utils_polynomial_selection.h"

void basic_polynomial_selection(polynomial_mpz * restrict polynomial, const mpz_t n, mpz_t m0, mpz_t m1, const unsigned long d)
{
    mpz_t d_root, tmp, tmp2, tmp3;
    mpz_inits(d_root, tmp, tmp2, tmp3, NULL);

    mpf_t tmpf;
    mpf_init(tmpf);

    mpf_set_prec(tmpf, 2048);

    mpf_set_z(tmpf, n);
    
    nth_root(tmpf, tmpf, d);

    mpz_set_f(d_root, tmpf);

    mpz_set(m0, d_root);
    mpz_set_ui(m1, 1);

    mpf_clear(tmpf);

    mpz_set(tmp, n);

    for (unsigned long i = 0 ; i <= d ; i++)
    {
        mpz_pow_ui(tmp2, d_root, d-i);

        mpz_div(tmp3, tmp, tmp2);

        set_coeff(polynomial, tmp3, d-i);

        mpz_mod(tmp, tmp, tmp2);
    }

    mpz_clears(d_root, tmp, tmp2, tmp3, NULL);
}

void Kleinjung_poly_selection(
    const mpz_t n,
    const dyn_array_classic * restrict primes,
    const unsigned long nb_roots,
    const unsigned long prime_bounds,
    const unsigned long multiplier,
    const mpz_t M,
    const unsigned long d,
    const unsigned long nb_poly_coarse_eval,
    const unsigned long nb_poly_precise_eval,
    gmp_randstate_t state,
    FILE *logfile
)
{
    dyn_array_classic kept_primes;
    init_classic(&kept_primes);

    for (size_t i = 0 ; i < primes->len ; i++)
    {
        if (primes->start[i] > prime_bounds) break;

        if (primes->start[i]%d == 1) append_classic(&kept_primes, primes->start[i]);
    }

    mpz_t a_d;
    mpz_init_set_ui(a_d, multiplier);

    mpf_t tmp_mpf, tmp_mpf2;
    mpf_inits(tmp_mpf, tmp_mpf2, NULL);

    mpf_set_prec(tmp_mpf, 2048);
    mpf_set_prec(tmp_mpf2, 2048);

    mpz_t tmp_mpz, a_d_max;
    mpz_inits(tmp_mpz, a_d_max, NULL);

    if (d < 4)
    {
        mpz_set(a_d_max, M);
    }
    else
    {
        mpz_pow_ui(tmp_mpz, M, 2*d - 2);
        mpf_set_z(tmp_mpf, tmp_mpz);
        mpf_set_z(tmp_mpf2, n);
        mpf_div(tmp_mpf, tmp_mpf, tmp_mpf2);

        nth_root(tmp_mpf2, tmp_mpf, d - 3);
        mpf_set_d(tmp_mpf, 0.5);
        mpf_add_ui(tmp_mpf2, tmp_mpf2, tmp_mpf);
        mpz_set_f(a_d_max, tmp_mpf2);
    }

    size_t cpt = 0;
    double avg = 0.0;

    mpz_t mw, ad1max, ad2max, tmp_mpz2;
    mpz_inits(mw, ad1max, ad2max, tmp_mpz2, NULL);

    while (mpz_cmp(a_d, a_d_max) < 0 && cpt < nb_poly_coarse_eval)
    {
        // Setup bounds

        mpf_set_z(tmp_mpf, n);
        mpf_set_z(tmp_mpf2, a_d);
        mpf_div(tmp_mpf, tmp_mpf, tmp_mpf2);
        nth_root(tmp_mpf2, tmp_mpf, d);

        mpz_set_f(mw, tmp_mpf2);
        mpz_add_ui(mw, mw, 1);


        mpz_mul(tmp_mpz, M, M);
        mpf_set_z(tmp_mpf, tmp_mpz);
        mpf_set_z(tmp_mpf2, mw);
        mpf_div(tmp_mpf, tmp_mpf, tmp_mpf2);

        mpf_set_d(tmp_mpf2, 0.5);
        mpf_add(tmp_mpf, tmp_mpf, tmp_mpf2);
        mpz_set_f(ad1max, tmp_mpf);


        if (d <= 2 || d == 4) mpz_set(ad2max, M);
        else if (d == 3) mpz_set(ad2max, mw);
        else
        {
            mpz_pow_ui(tmp_mpz, M, 2*d - 6);
            mpz_pow_ui(tmp_mpz2, mw, d - 4);

            mpf_set_z(tmp_mpf, tmp_mpz),
            mpf_set_z(tmp_mpf2, tmp_mpz2);

            mpf_div(tmp_mpf, tmp_mpf, tmp_mpf2);
            nth_root(tmp_mpf2, tmp_mpf, d - 2);

            mpz_set_f(ad2max, tmp_mpf2);
        }

        // Keep primes for which we have roots

        size_t nth_roots = 0;

        for (size_t i = 0 ; i < kept_primes.len ; i++)
        {
            unsigned long p = kept_primes.start[i];
            if (mpz_divisible_ui_p(a_d, p)) continue;

            mpz_set_ui(tmp_mpz, p);
            mpz_invert(tmp_mpz, a_d, tmp_mpz);
            mpz_mul(tmp_mpz, tmp_mpz, n);
            mpz_mod_ui(tmp_mpz, tmp_mpz, p);

            mpz_powm_ui(tmp_mpz2, tmp_mpz, (p - 1)/d, p);

            if (!mpz_cmp_ui(tmp_mpz2, 1)) nth_roots++;
        }

        unsigned long * Q = calloc(nth_roots, sizeof(unsigned long));
        unsigned long roots[nth_roots][d];

        unsigned long index = 0;

        polynomial_mpz tmp_poly;
        init_poly_degree(&tmp_poly, d);

        for (size_t i = 0 ; i < kept_primes.len ; i++)
        {
            unsigned long p = kept_primes.start[i];
            if (mpz_divisible_ui_p(a_d, p)) continue;

            set_coeff(&tmp_poly, a_d, d);

            mpz_neg(tmp_mpz, n);
            mpz_mod(tmp_mpz, tmp_mpz, p);

            set_coeff(&tmp_poly, tmp_mpz, 0);

            dyn_array_classic tmp_roots;
            init_classic(&tmp_roots);

            find_roots(&tmp_poly, &tmp_roots, p, state);

            if (tmp_roots.len)
            {
                Q[index] = p;
                for (size_t j = 0 ; j < d ; j++)
                {
                    roots[index][j] = tmo_roots.start[j];
                }

                index++;
            }

            free(tmp_roots.start);
        }

        free_polynomial(&tmp_poly);

        // Iterate over combinations of primes

        unsigned long * indexes = calloc(nb_roots, sizeof(unsigned long));
        for (size_t i = 0 ; i < nb_roots ; i++) indexes[i] = i;

        mpz_t product;
        mpz_init(product);

        while (true)
        {
            mpz_set_ui(product, 1);
            for (size_t i = 0 ; i < nb_roots ; i++) mpz_mul_ui(product, product, Q[indexes[i]]);

            if (mpz_cmp(product, ad1max) < 0)
            {
                // Do stuff
            }

            long j = nb_roots - 1;
            while (j >= 0 && indexes[j] == nth_roots - nb_roots + j)
                j--;

            if (j < 0)
                break;  // finished all combinations

            indexes[j]++;

            /* reset tail */
            for (size_t k = j + 1; k < nb_roots; k++) indexes[k] = indexes[k - 1] + 1;
        }

        free(indexes);
    }


    mpf_clears(tmp_mpf, tmp_mpf2, NULL);
    mpz_clears(tmp_mpz, tmp_mpz2, a_d, a_d_max, mw, ad1max, ad2max, product, NULL);
}