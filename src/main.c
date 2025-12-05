// Main file of the project

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gmp.h>

#include "logs.h"
#include "utils.h"
#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "NFS_relations.h"
#include "polynomial_functions.h"
#include "init_functions.h"
#include "generate_primes.h"
#include "polynomial_selection.h"
#include "algebraic_base.h"
#include "quadratic_characters.h"
#include "mono_cpu_sieve.h"
#include "build_matrix.h"
#include "block_lanczos.h"
#include "wiedemann.h"
#include "gaussian_elimination.h"
#include "extract_solution.h"

int main()
{
    printf("program started\n");

    // Setup logfile

    FILE *logfile = NULL;

    init_log(&logfile);

    // Initialize random seed

    srand(time(NULL));
    mpz_t random_a, random_b;
    mpz_init_set_ui(random_a, 3);
    mpz_init_set_ui(random_b, 82939);
    mpz_powm_ui(random_a, random_a, rand()*23%103, random_b);

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed(state, random_a);

    // Interactive input

    char number[100];
    printf("enter the number you wish to factor : ");
    fgets(number, 100, stdin);
    log_blank_line(logfile);

    mpz_t n;
    mpz_init(n);

    mpz_init_set_str(n, number, 10);

    // Initialization

    mpf_t ln10, ln2, e;
    mpf_inits(ln2, ln10, e, NULL);

    initialize_params(state, ln10, ln2, e);

    unsigned long degree = compute_degree(n, ln2, e);

    log_msg(logfile, "Degree selected is d = %lu", degree);
    log_blank_line(logfile);

    mpz_t smooth_bound;
    mpz_init(smooth_bound);

    compute_smooth_bound(n, smooth_bound, ln2, e);

    mpz_mul_ui(smooth_bound, smooth_bound, 4);

    // mpz_set_ui(smooth_bound, 256);

    mpz_t sieve_len;
    mpz_init_set(sieve_len, smooth_bound);
    mpz_mul_ui(sieve_len, sieve_len, 8);

    dyn_array_classic primes;
    init_classic(&primes);

    erasthotenes_sieve(&primes, smooth_bound);

    log_msg(logfile, "Rational Factor base of %lu primes generated.", primes.len);
    log_msg(logfile, "Largest prime = %lu.", primes.start[primes.len - 1]);

    dyn_array_classic logs;
    init_classic(&logs);

    compute_logs(&logs, primes);

    // Look for small factors

    for (size_t i = 0 ; i < primes.len ; i++)
    {
        if (mpz_divisible_ui_p(n, primes.start[i]))
        {
            mpz_t factor1, factor2;
            mpz_inits(factor1, factor2, NULL);

            char primality_factor1, primality_factor2;

            mpz_set_ui(factor1, primes.start[i]);
            mpz_divexact(factor2, n, factor1);

            if (mpz_probab_prime_p(factor1, 100) > 0)
            {
                primality_factor1 = 'p';
            } else {primality_factor1 = 'C';}

            if (mpz_probab_prime_p(factor2, 100) > 0)
            {
                primality_factor2 = 'p';
            } else {primality_factor2 = 'C';}

            log_blank_line(logfile);
            log_gmp_msg(logfile, "%Zd = %Zd (%c) x %Zd (%c)", n, factor1, primality_factor1, factor2, primality_factor2);
            if (logfile) fclose(logfile);

            mpz_clears(factor1, factor2);
            return 1;
        }
    }

    // The GNFS algorithm is split in 4 main steps


    // Polynomial selection
    // For now, very simple polynomial selection

    mpz_t m0, m1;
    mpz_inits(m0, m1, NULL);

    polynomial_mpz f_x;
    init_poly(&f_x);

    basic_polynomial_selection(&f_x, n, m0, m1, degree);

    mpz_t tmp;
    mpz_init(tmp);

    polynomial_mpz linear_poly;
    init_poly_degree(&linear_poly, 1);

    mpz_set(tmp, m1);
    set_coeff(&linear_poly, tmp, 1);

    mpz_set(tmp, m0);
    mpz_neg(tmp, tmp);
    set_coeff(&linear_poly, tmp, 0);

    mpz_t leading_coeff;
    mpz_init(leading_coeff);

    mpz_set(leading_coeff, f_x.coeffs[0]);

    polynomial_mpz g_x;
    init_poly_degree(&g_x, degree);

    mpz_set_ui(tmp, 1);
    set_coeff(&g_x, tmp, degree);

    set_coeff(&g_x, f_x.coeffs[1], degree - 1);

    for (size_t i = 2 ; i <= f_x.degree ; i++)
    {
        mpz_pow_ui(tmp, leading_coeff, i-1);
        mpz_mul(tmp, tmp, f_x.coeffs[i]);
        set_coeff(&g_x, tmp, f_x.degree - i);
    }

    mpz_clear(tmp);

    polynomial_mpz f_derivative;
    init_poly_degree(&f_derivative, degree - 1);
    poly_derivative(&f_derivative, f_x);

    polynomial_mpz g_derivative;
    init_poly_degree(&g_derivative, degree - 1);
    poly_derivative(&g_derivative, g_x);

    // Build algebraic factor base and quadratic characters

    // Algebraic factor base

    algebraic_base Algebraic_base;
    algebraic_base_init(&Algebraic_base);

    build_algebraic_base(&Algebraic_base, primes, f_x, n, state);

    size_t nb_Algebraic_pairs = 0;

    algebraic_base_prime *alg_prime = Algebraic_base.start;

    while (alg_prime != NULL) // Count the number of pairs in the Algebraic factor base
    {
        nb_Algebraic_pairs += alg_prime->roots.len;

        alg_prime = alg_prime->next;
    }

    // Large prime constants

    log_msg(logfile, "Algebraic base of size %zu generated.", nb_Algebraic_pairs);

    mpz_t large_prime_constant1, large_prime_constant2;
    mpz_inits(large_prime_constant1, large_prime_constant2, NULL);

    mpz_set_ui(large_prime_constant1, primes.start[primes.len - 1]);
    mpz_mul_ui(large_prime_constant1, large_prime_constant1, 100);

    mpz_mul_ui(large_prime_constant2, large_prime_constant1, primes.start[primes.len - 1]);

    // Quadratic characters base

    quadratic_character_base quad_char_base;
    quadratic_base_init(&quad_char_base);

    // unsigned long factor = create_quadratic_characters_base(&quad_char_base, f_x, f_derivative, n, leading_coeff, 1, mpz_get_ui(large_prime_constant2), state);
    unsigned long factor = create_quadratic_characters_base(&quad_char_base, f_x, f_derivative, n, leading_coeff, 3*mpz_sizeinbase(n, 2), mpz_get_ui(large_prime_constant2), state);

    if (factor) // If we have found a factor while building the quadratic characters base
    {
        mpz_t factor1, factor2;
        mpz_inits(factor1, factor2, NULL);

        char primality_factor1, primality_factor2;

        mpz_set_ui(factor1, factor);
        mpz_divexact(factor2, n, factor1);

        if (mpz_probab_prime_p(factor1, 100) > 0)
        {
            primality_factor1 = 'p';
        } else {primality_factor1 = 'C';}

        if (mpz_probab_prime_p(factor2, 100) > 0)
        {
            primality_factor2 = 'p';
        } else {primality_factor2 = 'C';}

        log_blank_line(logfile);
        log_gmp_msg(logfile, "%Zd = %Zd (%c) x %Zd (%c)", n, factor1, primality_factor1, factor2, primality_factor2);
        if (logfile) fclose(logfile);

        mpz_clears(factor1, factor2);
        return 1;
    }

    size_t nb_Quadratic_characters = 0;

    quadratic_character *quad_char = quad_char_base.start;

    while (quad_char != NULL) // Count the number of quadratic characters
    {
        nb_Quadratic_characters++;
        quad_char = quad_char->next;
    }

    log_msg(logfile, "Quadratic characters base of size %zu generated.", nb_Quadratic_characters);

    // Print and log stuff

    log_blank_line(logfile);
    log_gmp_msg(logfile, "Large prime bound 1 = %Zd = %lu*p_max", large_prime_constant1, primes.start[primes.len - 1]);
    log_gmp_msg(logfile, "Large prime bound 2 = %Zd = %lu*p_max^2", large_prime_constant2, primes.start[primes.len - 1]);
    log_blank_line(logfile);

    printf("f1(x) =\n");
    print_polynomial(&f_x);
    printf("\n");

    printf("f2(x) =\n");
    print_polynomial(&linear_poly);
    printf("\n");

    printf("g(x) =\n");
    print_polynomial(&g_x);
    printf("\n");

    // Compute some useful quantities

    polynomial_mpz tmp_poly, g_derivative_sq;
    init_poly(&tmp_poly);
    init_poly(&g_derivative_sq);

    poly_prod(&tmp_poly, g_derivative, g_derivative);
    poly_div(&g_derivative_sq, tmp_poly, g_x);

    free_polynomial(&tmp_poly);

    mpz_t g_derivative_eval;
    mpz_init(g_derivative_eval);

    mpz_mul(tmp, leading_coeff, m0);

    evaluate_homogeneous(g_derivative_eval, g_derivative, tmp, m1);
    mpz_mod(g_derivative_eval, g_derivative_eval, n);

    // Computing the set of small inert primes

    dyn_array_classic inert_set;
    init_classic(&inert_set);

    mpz_t tmp2;
    mpz_init(tmp2);

    for (size_t i = 0 ; i < primes.len ; i++)
    {
        mpz_mod_ui(tmp, m1, primes.start[i]);

        mpz_mod_ui(tmp2, leading_coeff, primes.start[i]);
        
        if (mpz_cmp_ui(tmp, 0) && mpz_cmp_ui(tmp2, 0) && irreducible(g_x, primes.start[i]))
        {
            append_classic(&inert_set, primes.start[i]);
        }
    }

    mpz_clear(tmp2);

    log_msg(logfile, "Found %lu inert primes in the factor base.", inert_set.len);
    log_blank_line(logfile);

    size_t len_divide_leading = 0;

    for (size_t i = 0 ; i < primes.len ; i++)
    {
        if (mpz_cmp_ui(leading_coeff, primes.start[i]) < 0)
        {
            break;
        }

        if (mpz_divisible_ui_p(leading_coeff, primes.start[i])) len_divide_leading++;
    }

    mpz_t *pow_div = calloc(len_divide_leading, sizeof(mpz_t));
    unsigned long *divide_leading = calloc(len_divide_leading, sizeof(unsigned long));

    size_t k = 0;

    for (size_t i = 0 ; i < primes.len ; i++)
    {
        if (mpz_cmp_ui(leading_coeff, primes.start[i]) < 0)
        {
            break;
        }

        if (mpz_divisible_ui_p(leading_coeff, primes.start[i]))
        {
            divide_leading[k] = primes.start[i];

            mpz_set_ui(tmp, primes.start[i]);
            while (mpz_divisible_p(leading_coeff, tmp)) mpz_mul_ui(tmp, tmp, primes.start[i]);

            mpz_divexact_ui(tmp, tmp, primes.start[i]);
            mpz_init_set(pow_div[k], tmp);

            k++;
        }
    }

    mpz_t prod_primes;
    mpz_init_set_ui(prod_primes, 1);
    for (size_t i = 0 ; i < primes.len ; i++) mpz_mul_ui(prod_primes, prod_primes, primes.start[i]);

    // Collecting relations

    nfs_relations relations;
    init_relations(&relations);

    mono_cpu_sieve(
        &relations,
        f_x,
        g_x,
        primes,
        Algebraic_base,
        nb_Algebraic_pairs,
        quad_char_base,
        nb_Quadratic_characters,
        leading_coeff,
        prod_primes,
        m0,
        m1,
        mpz_get_ui(sieve_len),
        large_prime_constant1,
        large_prime_constant2,
        divide_leading,
        pow_div,
        len_divide_leading,
        logs,
        state,
        logfile
    );

    // Linear algebra

    dyn_array_classic sparse_matrix;
    init_classic(&sparse_matrix);

    build_sparse_matrix(
        &sparse_matrix,
        &relations,
        &primes,
        &Algebraic_base,
        &quad_char_base,
        divide_leading,
        len_divide_leading
    );

    // Square root extraction

    bool *kernel_vector = calloc(relations.len, sizeof(bool));

    mpz_t divisor;
    mpz_init(divisor);

    while (true)
    {
        dyn_array kernel_vectors;
        init(&kernel_vectors);

        block_lanczos(&kernel_vectors, sparse_matrix, relations.len, 8, relations.len, logfile);

        for (size_t i = 0 ; i < kernel_vectors.len ; i++)
        {
            convert_to_vec(kernel_vectors.start[i], relations.len, kernel_vector);

            extract_solution(
                divisor,
                &relations,
                kernel_vector,
                &primes,
                g_x,
                g_derivative_sq,
                leading_coeff,
                n,
                m0,
                m1,
                g_derivative_eval,
                inert_set.start[inert_set.len-1],
                mpz_get_ui(sieve_len),
                state
            );

            mpz_gcd(divisor, n, divisor);

            mpz_divexact(tmp, n, divisor);

            gmp_printf("Factors found : %Zd %Zd\n\n", divisor, tmp);

            if (mpz_cmp_ui(divisor, 1) && mpz_cmp_ui(tmp, 1)) return 1;
        }

        free(kernel_vectors.start);
    }

    mpz_clears(n, m0, m1, divisor, NULL);
}