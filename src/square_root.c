#include <gmp.h>
#include <stdlib.h>

#include "dynamic_arrays.h"
#include "NFS_relations.h"
#include "polynomial_functions.h"
#include "square_root.h"

void extract_rational_square_root(
    mpz_t rational_square_root,
    nfs_relations *relations,
    bool *kernel_vector,
    mpz_t n,
    dyn_array_classic *rational_primes
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

    for (size_t i = 1 ; i < rational_primes->len ; i++)
    {
        mpz_set_ui(current_prime, rational_primes->start[i]);

        mpz_pow_ui(tmp_mpz, current_prime, prime_factors[i+1]>>1);

        mpz_mul(rational_square_root, rational_square_root, tmp_mpz);

        mpz_mod(rational_square_root, rational_square_root, n);
    }

    mpz_clears(tmp_mpz, current_prime, NULL);
}

void extract_algebraic_square_root(
    mpz_t algebraic_square_root,
    polynomial_mpz f_x,
    polynomial_mpz algebraic_square,
    mpz_t m0,
    mpz_t m1,
    mpz_t leading_coeff,
    mpz_t coeff_bound,
    unsigned long inert_prime,
    gmp_randstate_t state
)
{
    polynomial_mpz algebraic_root;
    init_poly(&algebraic_root);

    square_root_poly_mod(&algebraic_root, algebraic_square, f_x, inert_prime, state);

    printf("initial root computed\n");

    // Newton_lift();

    // evaluate_homogeneous(root);
}