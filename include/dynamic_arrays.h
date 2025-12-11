#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <gmp.h>
#include <stdbool.h>

typedef struct
{
    mpz_t* start;
    unsigned long len;
    unsigned long size;
    unsigned long initialized;
} dyn_array; // dynamic array of arbitrary large integers

typedef struct
{
    unsigned long *start;
    unsigned long len;
    unsigned long size;
} dyn_array_classic; // dynamic array of unsigned long integers

// Functions declaration

// Initialize functions

void init(dyn_array * restrict array);

void init_classic(dyn_array_classic * restrict array);

void init_len(dyn_array * restrict array, const unsigned long length);

void init_len_classic(dyn_array_classic * restrict array, const unsigned long length);

void init2_len(dyn_array * restrict array, const unsigned long length);

// Append functions

void append(dyn_array * restrict array, const mpz_t element);

void append_eco(dyn_array * restrict array, const mpz_t element);

void append_only(dyn_array * restrict array, const mpz_t element);

void append_only_si(dyn_array * restrict array, const signed long element);

void append_classic(dyn_array_classic * restrict array, const unsigned long element);

// Delete functions

void delete_classic(dyn_array_classic * restrict array, const unsigned long index);

void delete_dyn(dyn_array * restrict array, const unsigned long index);

void delete_dyn_eco(dyn_array * restrict array, const unsigned long index);

void delete_dyn_unsorted(dyn_array * restrict array, const unsigned long index);

// Insert functions

void insert_classic(dyn_array_classic * restrict array, const unsigned long element, const unsigned long index);

// Liberating arrays

void reset(dyn_array * restrict array);

void free_dyn_array(dyn_array * restrict array);

#endif // STRUCTURES_H