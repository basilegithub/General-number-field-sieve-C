#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

#include "dynamic_arrays.h"
#include "utils.h"

void init(dyn_array * restrict array)
{
	array->start = calloc(1, sizeof(mpz_t));
	array->len = 0;
	array->size = 1;
    array->initialized = 0;
}

void init_classic(dyn_array_classic * restrict array)
{
	array->start = calloc(1, sizeof(unsigned long));
	array->len = 0;
	array->size = 1;
}

void init_len(dyn_array * restrict array, const unsigned long length)
{
    array->size = 1;
	while (array->size <= length) array->size <<= 1;
	array->start = calloc(array->size, sizeof(mpz_t));
	array->len = length;
    array->initialized = 0;
}

void init_len_classic(dyn_array_classic * restrict array, const unsigned long length)
{
    array->size = 1;
	while (array->size <= length) array->size <<= 1;
	array->start = calloc(array->size, sizeof(unsigned long));
	array->len = length;
}

void init2_len(dyn_array * restrict array, const unsigned long length)
{
    array->size = 1;
	while (array->size <= length) array->size <<= 1;
	array->start = calloc(array->size, sizeof(mpz_t));
	array->len = length;
	for (unsigned long i = 0 ; i < length ; i++) mpz_init(array->start[i]);
    array->initialized = length;
}

void append(dyn_array * restrict array, const mpz_t element)
{
    if (array->size <= array->len)
    {
        array->size <<= 1;
        array->start = realloc(array->start, sizeof(mpz_t)*(array->size));
    }
    mpz_init_set(array->start[array->len], element);
    array->len++;
    array->initialized++;
}

void append_eco(dyn_array * restrict array, const mpz_t element)
{
    if (array->size <= array->len)
    {
        array->size <<= 1;
        array->start = realloc(array->start, sizeof(mpz_t)*(array->size));
    }

    if (array->len == array->initialized)
    {
        mpz_init_set(array->start[array->len], element);
        array->initialized++;
    }
    else
    {
        mpz_set(array->start[array->len], element);
    }
    array->len++;
}

void append_only(dyn_array * restrict array, const mpz_t element)
{
    mpz_set(array->start[array->len], element);
    array->len++;
}

void append_only_si(dyn_array * restrict array, const signed long element)
{
    mpz_set_si(array->start[array->len], element);
    array->len++;
}

void append_classic(dyn_array_classic * restrict array, const unsigned long element)
{
    if (array->size <= array->len)
    {
        array->size <<= 1;
        array->start = realloc(array->start, sizeof(unsigned long)*(array->size));
    }
    array->start[array->len] = element;
    array->len++;
}

void delete_classic(dyn_array_classic * restrict array, const unsigned long index)
{
    for (unsigned long i = index ; i < array->len-1 ; i++) array->start[i] = array->start[i+1];
    array->len--;
}

void delete_dyn(dyn_array * restrict array, const unsigned long index)
{
    for (unsigned long i = index ; i < array->len-1 ; i++) mpz_set(array->start[i], array->start[i+1]);
    array->len--;
    mpz_clear(array->start[array->len]);
    array->initialized--;
}

void delete_dyn_eco(dyn_array * restrict array, const unsigned long index)
{
    for (unsigned long i = index ; i < array->len-1 ; i++) mpz_set(array->start[i], array->start[i+1]);
    array->len--;
}

void delete_dyn_unsorted(dyn_array * restrict array, const unsigned long index)
{
    array->len--;
    mpz_set(array->start[index], array->start[array->len]);
    mpz_clear(array->start[array->len]);
    array->initialized--;
}

void insert_classic(dyn_array_classic * restrict array, const unsigned long element, const unsigned long index)
{
    if (array->size <= array->len)
    {
        array->size <<= 1;
        array->start = realloc(array->start, sizeof(unsigned long)*(array->size));
    }
    for (size_t i = array->len ; i > index ; i--)
    {
        array->start[i] = array->start[i-1];
    }
    array->start[index] = element;
    array->len++;
}

void reset(dyn_array * restrict array)
{
    array->len = 0;
}

void free_dyn_array(dyn_array * restrict array) {
    if (!array || !array->start) return;  // safety check

    // Clear each GMP integer
    for (size_t i = 0; i < array->initialized; i++) {
        mpz_clear(array->start[i]);
    }

    // Free the memory block holding the mpz_t structs
    free(array->start);

    // Reset the struct fields
    array->start = NULL;
    array->len = 0;
    array->size = 0;
    array->initialized = 0;
}