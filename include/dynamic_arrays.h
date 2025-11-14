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
    unsigned long* start;
    unsigned long len;
    unsigned long size;
} dyn_array_classic; // dynamic array of unsigned long integers

// Functions declaration

// Initialize functions

void init(dyn_array* array);
void init_classic(dyn_array_classic* array);
void init_len(dyn_array* array, unsigned long length);
void init_len_classic(dyn_array_classic* array, unsigned long length);
void init2_len(dyn_array* array, unsigned long length);

// Append functions

void append(dyn_array* array, mpz_t element);
void append_eco(dyn_array* array, mpz_t element);
void append_only(dyn_array* array, mpz_t element);
void append_only_si(dyn_array* array, signed long element);
void append_classic(dyn_array_classic* array, unsigned long element);

// Delete functions

void delete_classic(dyn_array_classic* array, unsigned long index);
void delete_classic_first(dyn_array_classic* array);
void delete_dyn(dyn_array* array, unsigned long index);
void delete_dyn_eco(dyn_array* array, unsigned long index);
void delete_dyn_unsorted(dyn_array* array, unsigned long index);

// Insert functions

void insert_classic(dyn_array_classic* array, unsigned long element, unsigned long index);

// Liberating arrays

void reset(dyn_array* array);
void free_dyn_array(dyn_array* array);

#endif // STRUCTURES_H