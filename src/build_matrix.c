#include <gmp.h>
#include <stdlib.h>

#include "dynamic_arrays.h"
#include "algebraic_base.h"
#include "quadratic_characters.h"
#include "NFS_relations.h"
#include "polynomial_functions.h"

void build_sparse_matrix(
    dyn_array_classic *sparse_matrix,
    nfs_relations *relations,
    dyn_array_classic *rational_primes,
    algebraic_base *algebraic_primes,
    quadratic_character_base *quad_char,
    unsigned long *divide_leading,
    unsigned long len_divide_leading
)
{
    mpz_t current_prime, tmp_mpz;
    mpz_inits(current_prime, tmp_mpz, NULL);

    unsigned long bit_count;
    bool parity_bit;


    // Encode rational side

    for (size_t i = 0 ; i < relations->len ; i++)
    {
        if (mpz_cmp_ui(relations->rels[i].rational_norm, 0) < 0) append_classic(sparse_matrix, i);
    }
    append_classic(sparse_matrix, relations->len);


    for (size_t i = 0 ; i < rational_primes->len ; i++)
    {
        mpz_set_ui(current_prime, rational_primes->start[i]);

        for (size_t j = 0 ; j < relations->len ; j++)
        {
            bit_count = mpz_remove(tmp_mpz, relations->rels[j].rational_norm, current_prime);
            parity_bit = bit_count&1;

            if (parity_bit)
            {
                append_classic(sparse_matrix, j);
            }
        }
        append_classic(sparse_matrix, relations->len);
    }


    // Encode algebraic side

    for (size_t i = 0 ; i < relations->len ; i++)
    {
        if (mpz_cmp_ui(relations->rels[i].algebraic_norm, 0) < 0) append_classic(sparse_matrix, i);
    }
    append_classic(sparse_matrix, relations->len);

    algebraic_base_prime *alg_prime = algebraic_primes->start;

    while (alg_prime != NULL)
    {
        mpz_set_ui(current_prime, alg_prime->prime);

        for (size_t i = 0 ; i < alg_prime->roots.len ; i++)
        {
            unsigned long root = alg_prime->roots.start[i];

            for (size_t j = 0 ; j < relations->len ; j++)
            {
                bit_count = mpz_remove(tmp_mpz, relations->rels[j].algebraic_norm, current_prime);
                parity_bit = bit_count&1;

                if (parity_bit)
                {
                    if (!evaluate_mod_p(relations->rels[j].poly_f, root, alg_prime->prime)) append_classic(sparse_matrix, j);
                }
            }

            append_classic(sparse_matrix, relations->len);
        }

        alg_prime = alg_prime->next;
    }

    // Encode quadratic characters side

    quadratic_character *quadratic = quad_char->start;

    while (quadratic != NULL)
    {
        unsigned long current_prime_ui = quadratic->q;
        mpz_set_ui(current_prime, quadratic->q);
        unsigned long root = quadratic->r;

        for (size_t i = 0 ; i < relations->len ; i++)
        {
            unsigned long res = evaluate_mod_p(relations->rels[i].poly_f, root, current_prime_ui);
            mpz_set_ui(tmp_mpz, res);

            if (mpz_legendre(tmp_mpz, current_prime) == -1)
            {
                append_classic(sparse_matrix, i);
            }
        }

        quadratic = quadratic->next;
    }

    // Encode b dividing the leading coefficient

    for (size_t i = 0 ; i < len_divide_leading ; i++)
    {
        mpz_set_ui(current_prime, divide_leading[i]);

        for (size_t j = 0 ; j < relations->len ; j++)
        {
            if (relations->rels[j].divide_leading[i])
            {
                bit_count = mpz_remove(tmp_mpz, relations->rels[j].algebraic_norm, current_prime);
                parity_bit = bit_count&1;

                if (parity_bit)
                {
                    append_classic(sparse_matrix, j);
                }
            }
            
            append_classic(sparse_matrix, relations->len);
        }
    }

    // Encode parity of used relations

    for (size_t i = 0 ; i < relations->len ; i++)
    {
        if (relations->rels[i].nb_relations&1) append_classic(sparse_matrix, i);
    }
    append_classic(sparse_matrix, relations->len);

    // Clear

    mpz_clears(current_prime, tmp_mpz, NULL);
}