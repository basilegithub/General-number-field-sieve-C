#include <gmp.h>

#include "dynamic_arrays.h"
#include "quadratic_characters.h"
#include "polynomial_functions.h"

void algebraic_base_init(quadratic_character_base *b)
{
    b->start = NULL;
    b->end = NULL;
}

void algebraic_base_clear(quadratic_character_base *b)
{
    quadratic_character *p = b->start;
    while (p != NULL) {
        quadratic_character *next = p->next;

        free(p);

        p = next;
    }

    b->start = NULL;
    b->end = NULL;
}

unsigned long create_quadratic_characters_base(quadratic_character_base *q_base, polynomial_mpz f, polynomial_mpz f_derivative, mpz_t n, mpz_t leading_coeff, unsigned long required_size, unsigned long start_prime)
{
    size_t cpt = 0;

    unsigned long test_number = start_prime + 2;
    if (!(test_number&1)) test_number--; // Ensure it starts at the first odd number after the last prime used in factor base, including large primes

    mpz_t tmp;
    mpz_init(tmp);

    while (cpt < required_size)
    {
        mpz_set_ui(tmp, test_number);
        if (mpz_probab_prime_p(tmp, 100) && mpz_divisible_ui_p(leading_coeff, test_number))
        {
            if (mpz_divisible_ui_p(n, test_number)) return test_number;

            dyn_array_classic roots;
            init_classic(&roots);

            basic_find_roots(f, &roots, test_number);
            for (size_t i = 0 ; i < roots.len ; i++)
            {
                if (!evaluate_mod_p(f_derivative, roots.start[i], test_number))
                {
                    quadratic_character *quad_char = malloc(sizeof(quadratic_character));
                    quad_char->q = test_number;
                    quad_char->r = roots.start[i];
                    quad_char->next = NULL;

                    if (q_base->start == NULL)
                    {
                        q_base->start = quad_char;
                    }
                    else
                    {
                        q_base->end->next = quad_char;
                    }

                    q_base->end = quad_char;

                    cpt++;
                }
            }

            free(roots.start);
        }

        test_number += 2;
    }

    mpz_clear(tmp);

    return 0;
}