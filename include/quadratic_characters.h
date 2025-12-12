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

void quadratic_base_init(quadratic_character_base * restrict b);

void quadratic_base_clear(quadratic_character_base * restrict b);

unsigned long create_quadratic_characters_base(
    quadratic_character_base * restrict q_base,
    const polynomial_mpz * restrict f,
    const polynomial_mpz * restrict f_derivative,
    const mpz_t n,
    const mpz_t leading_coeff,
    const unsigned long required_size,
    const unsigned long start_prime, 
    gmp_randstate_t state
);

#endif // QUADRATIC_CHARACTERS_H