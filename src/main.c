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
#include "init_functions.h"
#include "generate_primes.h"
#include "polynomial_selection.h"

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

    mpz_t n;
    mpz_init(n);

    mpz_init_set_str(n, number, 10);

    // Initialization

    mpf_t ln10, ln2, e;
    mpf_inits(ln2, ln10, e, NULL);

    initialize_params(state, ln10, ln2, e);

    unsigned long degree = compute_degree(n, ln2, e);

    log_msg(logfile, "Degree selected is d = %lu", degree);

    mpz_t smooth_bound;
    mpz_init(smooth_bound);

    compute_smooth_bound(n, smooth_bound, ln2, e);

    mpz_set_ui(smooth_bound, 10);

    dyn_array_classic primes;
    init_classic(&primes);

    erasthotenes_sieve(&primes, smooth_bound);

    log_msg(logfile, "Factor base of %lu primes generated.", primes.len);
    log_msg(logfile, "Largest prime = %lu", primes.start[primes.len - 1]);

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

    mpz_t large_prime_constant1, large_prime_constant2;
    mpz_inits(large_prime_constant1, large_prime_constant2, NULL);

    mpz_set_ui(large_prime_constant1, primes.start[primes.len - 1]);
    mpz_mul_ui(large_prime_constant1, large_prime_constant1, 100);

    mpz_mul_ui(large_prime_constant2, large_prime_constant1, primes.start[primes.len - 1]);

    log_gmp_msg(logfile, "Large prime bound 1 = %Zd = %lu*p_max", large_prime_constant1, primes.start[primes.len - 1]);
    log_gmp_msg(logfile, "Large prime bound 2 = %Zd = %lu*p_max^2", large_prime_constant2, primes.start[primes.len - 1]);

    polynomial_mpz polynomial;
    init_poly(&polynomial);

    basic_polynomial_selection(&polynomial, n, degree);

    print_polynomial(&polynomial);

    // The GNFS algorithm is split in 4 main steps

    // Polynomial selection

    // For now, very easy polynomial selection

    // Sieving

    // Linear algebra

    // Square root extraction



    mpz_clear(n);
}