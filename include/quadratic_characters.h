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
} quadratic_character_base; // Dynamic algebraic base

// Functions definition

void quadratic_character_base_init(quadratic_character_base *b);
void algebraic_base_clear(algebraic_base *b);
void build_algebraic_base(algebraic_base *b, dyn_array_classic primes, polynomial_mpz g_x, mpz_t n);

#endif // QUADRATIC_CHARACTERS_H