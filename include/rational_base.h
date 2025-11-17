#ifndef RATIONAL_BASE_H
#define RATIONAL_BASE_H

#include <gmp.h>

#include "dynamic_arrays.h"

typedef struct
{
    unsigned long prime;
    dyn_array_classic roots;
    struct rational_base_prime *next;
} rational_base_prime; // Dynamic rational base prime

typedef struct
{
    rational_base_prime *start;
    rational_base_prime *end;
} rational_base; // Dynamic rational base

// Functions definition

void rational_base_init(rational_base *b);
void rational_base_clear(rational_base *b);

#endif // RATIONAL_BASE_H