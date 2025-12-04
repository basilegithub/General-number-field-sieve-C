#include <gmp.h>
#include <stdlib.h>

#include "dynamic_arrays.h"
#include "NFS_relations.h"

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