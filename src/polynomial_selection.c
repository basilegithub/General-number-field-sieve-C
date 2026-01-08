#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "utils.h"
#include "utils_polynomial_selection.h"
#include "polynomial_selection.h"

void init_ranking(
    polynomial_ranking * restrict ranking,
    size_t size
)
{
    ranking->ranking = calloc(size, sizeof(polynomial_ranking_element));
    ranking->size = size;
    ranking->len = 0;
}

void append_element(
    polynomial_ranking * restrict ranking,
    polynomial_ranking_element * restrict element
)
{
    // It is assumed that we have checked that the ranking still have space to append an element
    init_poly(&ranking->ranking[ranking->len].poly);
    copy_polynomial(&ranking->ranking[ranking->len].poly, &element->poly);

    mpz_init_set(ranking->ranking[ranking->len].m0, element->m0);
    mpz_init_set(ranking->ranking[ranking->len].m1, element->m1);

    ranking->len++;
}

void insert_element(
    polynomial_ranking * restrict ranking,
    polynomial_ranking_element * restrict element,
    size_t index
)
{
    // If there is enough space, the ranking is extended.
    // If not, the last polynomial is discarded

    if (ranking->len < ranking->size)
    {
        append_element(ranking, ranking->ranking[ranking->len - 1]); // Shift last element

        for (size_t i = ranking->len - 2 ; i > index ; i--)
        {
            copy_polynomial(&ranking->ranking[i].poly, &ranking->ranking[i - 1].poly);
            mpz_set(ranking->ranking[i].m0, ranking->ranking[i - 1].m0);
            mpz_set(ranking->ranking[i].m1, ranking->ranking[i - 1].m1);
        }
    }
    else
    {
        for (size_t i = ranking->len - 1 ; i > index ; i--) // Element at ranking->len - 1 will be overriden
        {
            copy_polynomial(&ranking->ranking[i].poly, &ranking->ranking[i - 1].poly);
            mpz_set(ranking->ranking[i].m0, ranking->ranking[i - 1].m0);
            mpz_set(ranking->ranking[i].m1, ranking->ranking[i - 1].m1);
        }
    }

    copy_polynomial(&ranking->ranking[index].poly, &element->poly);

    mpz_set(ranking->ranking[index].m0, element->m0);
    mpz_set(ranking->ranking[index].m1, element->m1);

    ranking->len++;
}

void free_ranking(
    polynomial_ranking * restrict ranking
)
{
    for (size_t i = 0 ; i < ranking->len ; i++)
    {
        free_polynomial(&ranking->ranking[i].poly);
        mpz_clear(ranking->ranking[i].m0);
        mpz_clear(ranking->ranking[i].m1);
    }

    free(ranking->ranking);
    ranking->len = 0;
    ranking->size = 0;
}

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
    polynomial_mpz * restrict res_poly,
    mpz_t m0,
    mpz_t m1,
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
        mpf_add(tmp_mpf2, tmp_mpf2, tmp_mpf);
        mpz_set_f(a_d_max, tmp_mpf2);
    }

    size_t cpt = 0;
    double avg = 0.0;

    mpz_t mw, ad1max, ad2max, tmp_mpz2;
    mpz_inits(mw, ad1max, ad2max, tmp_mpz2, NULL);

    polynomial_mpz tmp_poly;
    init_poly(&tmp_poly);

    unsigned long * tmp_indexes = calloc(nb_roots, sizeof(unsigned long));

    mpz_t tmp_m_mu, tmp_prod;
    mpz_inits(tmp_m_mu, tmp_prod, NULL);

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

            mpz_set_ui(tmp_mpz2, p);

            mpz_powm_ui(tmp_mpz2, tmp_mpz, (p - 1)/d, tmp_mpz2);

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
            mpz_mod_ui(tmp_mpz, tmp_mpz, p);

            set_coeff(&tmp_poly, tmp_mpz, 0);

            dyn_array_classic tmp_roots;
            init_classic(&tmp_roots);

            find_roots(&tmp_poly, &tmp_roots, p, state);

            if (tmp_roots.len)
            {
                Q[index] = p;
                for (size_t j = 0 ; j < d ; j++)
                {
                    roots[index][j] = tmp_roots.start[j];
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

        mpz_t roots_used[nb_roots][d];
        for (size_t i = 0 ; i < nb_roots ; i++)
        {
            for (size_t j = 0 ; j < d ; j++)
            {
                mpz_init(roots_used[i][j]);
            }
        }

        mpz_t m0_local;
        mpz_init(m0_local);

        while (true)
        {
            mpz_set_ui(product, 1);
            for (size_t i = 0 ; i < nb_roots ; i++) mpz_mul_ui(product, product, Q[indexes[i]]);

            if (mpz_cmp(product, ad1max) < 0)
            {
                dyn_array_classic Q_used;
                init_classic(&Q_used);

                for (size_t i = 0 ; i < nb_roots ; i++)
                {
                    append_classic(&Q_used, Q[indexes[i]]);

                    for (size_t j = 0 ; j < d ; j++)
                    {
                        mpz_set_ui(roots_used[i][j], roots[indexes[i]][j]);
                    }
                }

                for (size_t i = 0 ; i < nb_roots ; i++)
                {
                    mpz_divexact_ui(tmp_mpz, product, Q_used.start[i]);
                    mpz_set_ui(tmp_mpz2, Q_used.start[i]);
                    mpz_invert(tmp_mpz2, tmp_mpz, tmp_mpz2);
                    mpz_mul(tmp_mpz, tmp_mpz, tmp_mpz2);

                    for (size_t j = 0 ; j < d ; j++)
                    {
                        mpz_mul(roots_used[i][j], roots_used[i][j], tmp_mpz);
                        mpz_mod(roots_used[i][j], roots_used[i][j], product);
                    }
                }

                mpz_neg(m0_local, mw);
                mpz_mod(m0_local, m0_local, product);
                mpz_sub(m0_local, mw, m0_local);

                mpz_t e_array[nb_roots][d];

                for (size_t i = 0 ; i < nb_roots ; i++)
                {
                    for (size_t j = 0 ; j < d ; j++)
                    {
                        mpz_init(e_array[i][j]);
                    }
                }

                compute_e_array(
                    m0_local,
                    nb_roots,
                    d,
                    roots_used,
                    product,
                    a_d,
                    n,
                    e_array
                );

                mpf_t f0;
                mpf_init(f0);

                mpf_t * f = calloc(nb_roots*d, sizeof(mpf_t));

                for (size_t i = 0 ; i < nb_roots*d ; i++)
                {
                    mpf_init(f[i]);
                }

                compute_f(
                    n,
                    a_d,
                    m0_local,
                    d,
                    product,
                    nb_roots,
                    roots_used,
                    e_array,
                    f,
                    f0
                );

                mpf_t epsilon;
                mpf_init(epsilon);

                mpf_set_z(tmp_mpf, ad2max);
                mpf_set_z(tmp_mpf2, m0_local);
                mpf_div(epsilon, tmp_mpf, tmp_mpf2);

                unsigned long len_vec1 = nb_roots>>1;

                unsigned long nrows1 = my_power(d, len_vec1);
                unsigned long ncols1 = len_vec1;

                unsigned long * array1_indices = calloc(nrows1*ncols1, sizeof(unsigned long));

                mpf_t * array1_u_values = calloc(nrows1, sizeof(mpf_t));

                for (size_t i = 0 ; i < nrows1 ; i++)
                {
                    mpf_init(array1_u_values[i]);
                }

                create_first_array(
                    nb_roots,
                    d,
                    f0,
                    f,
                    nrows1,
                    ncols1,
                    array1_u_values,
                    array1_indices
                );

                unsigned long len_vec2 = nb_roots - len_vec1;

                unsigned long nrows2 = my_power(d, len_vec2);
                unsigned long ncols2 = len_vec2;

                unsigned long * array2_indices = calloc(nrows2*ncols2, sizeof(unsigned long));

                mpf_t * array2_u_values = calloc(nrows2, sizeof(mpf_t));

                for (size_t i = 0 ; i < nrows2 ; i++)
                {
                    mpf_init(array2_u_values[i]);
                }

                create_second_array(
                    nb_roots,
                    len_vec1,
                    d,
                    nrows2,
                    ncols2,
                    f,
                    array2_u_values,
                    array2_indices
                );

                size_t minimum = 0;

                for (size_t i = 0 ; i < nrows2 ; i++)
                {
                    mpf_add(tmp_mpf, array1_u_values[minimum], epsilon);
                    while (minimum + 1 < nrows1 && mpf_cmp(array2_u_values[i], tmp_mpf) > 0)
                    {
                        minimum++;
                        mpf_add(tmp_mpf, array1_u_values[minimum], epsilon);
                    }

                    if (minimum == nrows1) break;

                    size_t z = minimum;
                    if (z < nrows1)
                    {
                        mpf_sub(tmp_mpf, array2_u_values[i], array1_u_values[z]);
                        mpf_abs(tmp_mpf, tmp_mpf);
                        while (z < nrows1 && mpf_cmp(tmp_mpf, epsilon) < 0)
                        {
                            // Compute the polynomial, the associated m_mu and prod
                            for (size_t i = 0 ; i < len_vec1 ; i++) tmp_indexes[i] = array1_indices[i];
                            for (size_t i = 0 ; i < len_vec2 ; i++) tmp_indexes[len_vec1 + i] = array2_indices[i];
                            mpz_set(tmp_prod, product);

                            compute_m_mu(tmp_m_mu, m0_local, nb_roots, d, roots_used, tmp_indexes);
                            build_poly_coeffs(&tmp_poly, tmp_m_mu, tmp_prod, a_d, n, d);
                            // Get the L2 norm value of the generated polynomial
                            // Rank the polynomial, keep it if good
                            // If enough polynomial have been checked, finish


                            z++;
                        }
                    }

                }



                for (size_t i = 0 ; i < nb_roots ; i++)
                {
                    for (size_t j = 0 ; j < d ; j++)
                    {
                        mpz_clear(e_array[i][j]);
                    }
                }

                for (size_t i = 0 ; i < nb_roots*d ; i++)
                {
                        mpf_clear(f[i]);
                }

                free(f);

                for (size_t i = 0 ; i < nrows1 ; i++)
                {
                    mpf_clear(array1_u_values[i]);
                }

                for (size_t i = 0 ; i < nrows2 ; i++)
                {
                    mpf_clear(array2_u_values[i]);
                }

                free(array1_u_values);
                free(array2_u_values);
                free(array1_indices);
                free(array2_indices);

                free(Q_used.start);

                mpf_clears(f0, epsilon, NULL);
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

        free(Q);

        mpz_clears(product, m0_local, NULL);
    }

    free(kept_primes.start);
    free(tmp_indexes);

    free_polynomial(&tmp_poly);

    mpf_clears(tmp_mpf, tmp_mpf2, NULL);
    mpz_clears(tmp_mpz, tmp_mpz2, a_d, a_d_max, mw, ad1max, ad2max, tmp_m_mu, tmp_prod, NULL);
}