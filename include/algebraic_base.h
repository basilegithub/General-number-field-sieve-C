#ifndef ALGEBRAIC_BASE_H_BASE_H
#define ALGEBRAIC_BASE_H_BASE_H

#include <gmp.h>

#include "dynamic_arrays.h"

typedef struct
{
    unsigned long prime;
    dyn_array_classic roots;
    struct algebraic_base_prime *next;
} algebraic_base_prime; // Dynamic algebraic base prime

typedef struct
{
    algebraic_base_prime *start;
    algebraic_base_prime *end;
} algebraic_base; // Dynamic algebraic base

// Functions definition

void algebraic_base_init(algebraic_base *b);
void algebraic_base_clear(algebraic_base *b);

#endif // ALGEBRAIC_BASE_H