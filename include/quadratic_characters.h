#ifndef QUADRATIC_CHARACTERS_H
#define QUADRATIC_CHARACTERS_H

#include <gmp.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"

typedef struct quadratic_character
{
    unsigned long q;
    unsigned long r;
    struct quadratic_character *next;
} quadratic_character; // Dynamic algebraic base prime

typedef struct
{
    quadratic_character *start;
    quadratic_character *end;
} quadratic_character_base; // Dynamic algebraic base

// Functions definition

void quadratic_character_base_init(quadratic_character_base *b);
void algebraic_base_clear(quadratic_character_base *b);
void create_quadratic_characters_base(quadratic_character_base *q_base, polynomial_mpz f, polynomial_mpz f_derivative, mpz_t n, mpz_t leading_coeff, unsigned long required_size, unsigned long start_prime);

#endif // QUADRATIC_CHARACTERS_H