#ifndef SMOOTH_TEST_H
#define SMOOTH_TEST_H

#include <gmp.h>
#include <stdbool.h>

#include "dynamic_arrays.h"
#include "single_linked_list.h"
#include "polynomial_structures.h"
#include "NFS_relations.h"
#include "utils.h"


void pollard_rho(const mpz_t m, mpz_t p1, mpz_t p2, gmp_randstate_t state)
{
    mpz_t a, b, d, e, x, y, X, Y, tmp;
    mpz_inits(a, b, d, e, x, y, X, Y, tmp, NULL);

    while (true)
    {
        mpz_sub_ui(tmp, m, 4);
        mpz_urandomm(a, state, tmp);
        mpz_add_ui(a, a, 1);

        mpz_sub_ui(tmp, m, 1);
        mpz_urandomm(b, state, tmp);

        mpz_set(x, b);
        mpz_set(y, b);
        mpz_set_ui(d, 1);

        while (mpz_cmp_ui(d, 1) == 0)
        {
            mpz_set_ui(e, 1);
            mpz_set(X, x);
            mpz_set(Y, y);

            for (unsigned int k = 0 ; k < 100 ; k++)
            {
                mpz_mul(tmp, x, x);
                mpz_add(tmp, tmp, a);
                mpz_mod(x, tmp, m); // x = (x*x+a)%m

                mpz_mul(tmp, y, y);
                mpz_add(tmp, tmp, a);
                mpz_mod(y, tmp, m); // y = (y*y+a)%m

                mpz_mul(tmp, y, y);
                mpz_add(tmp, tmp, a);
                mpz_mod(y, tmp, m); // y = (y*y+a)%m

                mpz_sub(tmp, x, y);
                mpz_mod(tmp, tmp, m);
                mpz_mul(b, e, tmp);
                mpz_mod(e, b, m); // e = e*(x-y)%m
            }

            mpz_gcd(d, e, m);
        }

        if (mpz_cmp(d, m) == 0)
        {
            mpz_set(x, X);
            mpz_set(y, Y);
            mpz_set_ui(d, 1);

            while (mpz_cmp_ui(d, 1) == 0)
            {
                mpz_mul(tmp, x, x);
                mpz_add(tmp, tmp, a);
                mpz_mod(x, tmp, m); // x = (x*x+a)%m

                mpz_mul(tmp, y, y);
                mpz_add(tmp, tmp, a);
                mpz_mod(y, tmp, m); // y = (y*y+a)%m

                mpz_mul(tmp, y, y);
                mpz_add(tmp, tmp, a);
                mpz_mod(y, tmp, m); // y = (y*y+a)%m

                mpz_sub(tmp, x, y);
                mpz_mod(tmp, tmp, m);
                mpz_gcd(d, tmp, m);
            }
        }

        if (mpz_cmp(d, m)) // If d != m
        {
            mpz_divexact(tmp, m, d);
            if (mpz_cmp(tmp, d) < 0)
            {
                mpz_set(p1, tmp);
                mpz_set(p2, d);
            }
            else
            {
                mpz_set(p1, d);
                mpz_set(p2, tmp);
            }
            mpz_clears(a, b, d, e, x, y, X, Y, tmp, NULL);
            return;
        }
    }
}

unsigned long build_product_tree(
    const dyn_array * restrict reported,
    dyn_array * restrict tree_array,
    const mpz_t prod_primes,
    const mpz_t prod_primes_p1,
    mpz_t tmp
)
{
    unsigned long tmp_long = reported->len;
    unsigned long index1 = tmp_long, index2 = tmp_long;

    tmp_long >>= 1;
    for (unsigned long i = 0 ; i < tmp_long ; i++)
    {
        mpz_mul(tmp, reported->start[i<<1], reported->start[(i<<1)+1]);
        if (mpz_cmp(tmp, prod_primes) < 1) mpz_set(tree_array->start[index1+i], tmp);
        else mpz_set(tree_array->start[index1+i], prod_primes_p1);
    }

    index2 += tmp_long;
    
    while (tmp_long > 1)
    {
        for (unsigned long i = 0 ; i < tmp_long ; i++)
        {
            if (!mpz_cmp(tree_array->start[index1+(i<<1)], prod_primes_p1) || !mpz_cmp(tree_array->start[index1+(i<<1)+1], prod_primes_p1))
            {
                mpz_set(tree_array->start[index2+i], prod_primes_p1);
            }
            else
            {
                mpz_mul(tmp, tree_array->start[index1+(i<<1)], tree_array->start[index1+(i<<1)+1]);

                if (mpz_cmp(tmp, prod_primes) < 1) mpz_set(tree_array->start[index2+i], tmp);
                else mpz_set(tree_array->start[index2+i], prod_primes_p1);
            }
        }
        index1 = index2;
        tmp_long >>= 1;
        index2 += tmp_long;
    }
    return tmp_long;
}

void build_remainder_tree(
    const dyn_array * restrict reported,
    dyn_array * restrict tree_array,
    const mpz_t prod_primes,
    const mpz_t prod_primes_p1,
    unsigned long tmp_long
)
{
    unsigned long index1 = 0;
    unsigned long index2 = 1;
    unsigned long stored_len = tree_array->len;
    mpz_mod(tree_array->start[stored_len-1], prod_primes, tree_array->start[stored_len-1]);

    bool flag_overflow;

    while (tmp_long < reported->len>>1)
    {
        for (unsigned long i = 0 ; i < tmp_long ; i++)
        {
            flag_overflow = !mpz_cmp(tree_array->start[stored_len-1-index2-(i<<1)], prod_primes_p1);
            if (mpz_cmp(tree_array->start[stored_len-1-index2-(i<<1)], prod_primes_p1))
            {
                mpz_mod(tree_array->start[stored_len-1-index2-(i<<1)], tree_array->start[stored_len-1-index1-i], tree_array->start[stored_len-1-index2-(i<<1)]);
            }
            else mpz_set(tree_array->start[stored_len-1-index2-(i<<1)], prod_primes);

            if (flag_overflow)
            {
                mpz_mod(tree_array->start[stored_len-2-index2-(i<<1)], tree_array->start[stored_len-1-index1-i], tree_array->start[stored_len-2-index2-(i<<1)]);
            }
            else mpz_set(tree_array->start[stored_len-2-index2-(i<<1)], prod_primes);
        }
        index1 = index2;
        tmp_long <<= 1;
        index2 += tmp_long;
    }


    for (unsigned long i = 0 ; i < tmp_long ; i++)
    {
        mpz_mod(tree_array->start[i<<1], tree_array->start[reported->len+i], reported->start[i<<1]);
        mpz_mod(tree_array->start[(i<<1)+1], tree_array->start[reported->len+i], reported->start[(i<<1)+1]);
    }
}

void batch_smooth(
    nfs_relations * restrict smooth_candidates,
    dyn_array *reported,
    dyn_array *tree_array,
    const mpz_t prod_primes,
    const mpz_t prod_primes_p1,
    const mpz_t limit,
    const mpz_t limit_2,
    const size_t index,
    gmp_randstate_t state
)
{
    mpz_t tmp, tmp2;
    mpz_inits(tmp, tmp2, NULL);
    mpz_set_ui(tmp, 2);

    for (size_t i = 0 ; i < reported->len ; i++)
    {
        mpz_mul(reported->start[i], smooth_candidates->rels[index + i].rational_norm, smooth_candidates->rels[index + i].algebraic_norm);
        mpz_abs(reported->start[i], reported->start[i]);
        mpz_remove(reported->start[i], reported->start[i], tmp);
    } // Remove occurences of 2

    unsigned long squaring_steps = 0;
    
    for (mpz_t* ptr = reported->start ; ptr < reported->start + reported->len ; ptr++)
    {
        while (mpz_cmp(tmp, *ptr) < 0)
        {
            mpz_mul(tmp2, tmp, tmp);
            mpz_set(tmp, tmp2);
            squaring_steps++;
        }
    }
    // computed log(log(max[reported]))

    // computing prod tree
    unsigned long tmp_long = build_product_tree(reported, tree_array, prod_primes, prod_primes_p1, tmp);

    //computing remainder tree
    build_remainder_tree(reported, tree_array, prod_primes, prod_primes_p1, tmp_long);

    mpz_t remaining_factor, prime_divisor_1, prime_divisor_2;
    mpz_inits(remaining_factor, prime_divisor_1, prime_divisor_2, NULL);

    for (unsigned long i = 0 ; i < reported->len ; i++)
    {
        tmp_long = 0;
        mpz_set(tmp, tree_array->start[i]);
        while (mpz_cmp_ui(tmp, 0) && tmp_long < squaring_steps)
        {
            mpz_mul(tmp, tmp, tmp);
            mpz_mod(tmp, tmp, reported->start[i]);
            tmp_long++; 
        }

        mpz_set(remaining_factor, tmp);

        if (!mpz_cmp_ui(remaining_factor, 0))
        {
            continue;
        }
        else
        {
            mpz_gcd(tmp2, reported->start[i], remaining_factor);
            mpz_divexact(remaining_factor, reported->start[i], tmp2);

            if (mpz_cmp(remaining_factor, limit) < 1)
            {
                if (mpz_divisible_p(smooth_candidates->rels[index + i].rational_norm, tmp2))
                {
                    append_classic(&smooth_candidates->rels[index + i].rational_large_primes, mpz_get_ui(tmp2));
                }
                else
                {
                    list_append(&smooth_candidates->rels[index + i].algebraic_large_primes, mpz_get_ui(tmp2), 1);
                }
            }
            else if (mpz_cmp(remaining_factor, limit_2) < 1 && !fermat_primality(remaining_factor))
            {
            pollard_rho(remaining_factor, prime_divisor_1, prime_divisor_2, state);

            if (mpz_divisible_p(smooth_candidates->rels[index + i].rational_norm, remaining_factor))
            {
                append_classic(&smooth_candidates->rels[index + i].rational_large_primes, mpz_get_ui(prime_divisor_1));
                append_classic(&smooth_candidates->rels[index + i].rational_large_primes, mpz_get_ui(prime_divisor_2));
            }
            else
            {
                list_append(&smooth_candidates->rels[index + i].algebraic_large_primes, mpz_get_ui(prime_divisor_1), 1);
                list_append(&smooth_candidates->rels[index + i].algebraic_large_primes, mpz_get_ui(prime_divisor_2), 1); 
            }
            }
            else smooth_candidates->rels[index + i].nb_relations = 0;
        }
    }
    mpz_clears(remaining_factor, prime_divisor_1, prime_divisor_2, tmp, tmp2, NULL);
}

void naive_smooth(
    nfs_relations * restrict smooth_candidates,
    const dyn_array_classic * restrict primes,
    const mpz_t limit,
    const mpz_t limit_2,
    const size_t index,
    gmp_randstate_t state
)
{
    mpz_t remaining_factor, prime_divisor_1, prime_divisor_2;
    mpz_inits(remaining_factor, prime_divisor_1, prime_divisor_2, NULL);

    unsigned long div_count;
    const unsigned long max_prime = primes->start[primes->len-1];
    bool flag_is_smooth;

    for (size_t i = index ; i < smooth_candidates->len ; i++)
    {
        mpz_mul(remaining_factor, smooth_candidates->rels[i].rational_norm, smooth_candidates->rels[i].algebraic_norm);
        mpz_abs(remaining_factor, remaining_factor);
        flag_is_smooth = false;
        
        for (unsigned long j = 0 ; j < primes->len ; j++)
        {
            mpz_set_ui(prime_divisor_1, primes->start[j]);
            div_count = mpz_remove(remaining_factor, remaining_factor, prime_divisor_1);

            if (div_count > 0 && mpz_cmp_ui(remaining_factor, max_prime) < 1)
            {
                flag_is_smooth = true;
                break;
            }
        }

        if (flag_is_smooth)
        {
            continue; // Do nothing
        }
        else if (mpz_cmp(remaining_factor, limit) < 1)
        {
            if (mpz_divisible_p(smooth_candidates->rels[i].rational_norm, remaining_factor))
            {
                append_classic(&smooth_candidates->rels[i].rational_large_primes, mpz_get_ui(remaining_factor));
            }
            else
            {
                list_append(&smooth_candidates->rels[i].algebraic_large_primes, mpz_get_ui(remaining_factor), 1);
            }
        }
        else if (mpz_cmp(remaining_factor, limit_2) < 1 && !fermat_primality(remaining_factor))
        {
            pollard_rho(remaining_factor, prime_divisor_1, prime_divisor_2, state);

            if (mpz_divisible_p(smooth_candidates->rels[i].rational_norm, remaining_factor))
            {
                append_classic(&smooth_candidates->rels[i].rational_large_primes, mpz_get_ui(prime_divisor_1));
                append_classic(&smooth_candidates->rels[i].rational_large_primes, mpz_get_ui(prime_divisor_2));
            }
            else
            {
                list_append(&smooth_candidates->rels[i].algebraic_large_primes, mpz_get_ui(prime_divisor_1), 1);
                list_append(&smooth_candidates->rels[i].algebraic_large_primes, mpz_get_ui(prime_divisor_2), 1); 
            }
        }
        else smooth_candidates->rels[i].nb_relations = 0; // Not smooth

    }
    mpz_clears(remaining_factor, prime_divisor_1, prime_divisor_2, NULL);
}


#endif // SMOOTH_TEST_H