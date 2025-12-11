#include <stdlib.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "algebraic_base.h"
#include "polynomial_functions.h"

void algebraic_base_init(algebraic_base *b)
{
    b->start = NULL;
    b->end = NULL;
}

void algebraic_base_clear(algebraic_base *b)
{
    algebraic_base_prime *p = b->start;
    while (p != NULL) {
        algebraic_base_prime *next = p->next;

        free(&p->roots);

        free(p);

        p = next;
    }

    b->start = NULL;
    b->end = NULL;
}

void build_algebraic_base(algebraic_base *b, const dyn_array_classic * restrict primes, const polynomial_mpz g_x, gmp_randstate_t state)
{
    algebraic_base_prime *alg_prime = malloc(sizeof(algebraic_base_prime));
    init_classic(&alg_prime->roots);

    unsigned long first_prime = primes->start[0];

    // basic_find_roots(g_x, &alg_prime->roots, first_prime);
    find_roots(g_x, &alg_prime->roots, first_prime, state);

    alg_prime->prime = first_prime;
    alg_prime->next = NULL;

    b->start = alg_prime;
    b->end = alg_prime;

    for (size_t i = 1 ; i < primes->len ; i++)
    {
        algebraic_base_prime *next = malloc(sizeof(algebraic_base_prime));
        init_classic(&next->roots);

        unsigned long next_prime = primes->start[i];

        // basic_find_roots(g_x, &next->roots, next_prime);
        find_roots(g_x, &next->roots, next_prime, state);

        next->prime = next_prime;
        next->next = NULL;

        b->end->next = next;
        b->end = next;
    }
}