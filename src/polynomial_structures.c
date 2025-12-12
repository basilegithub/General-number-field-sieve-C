#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "polynomial_structures.h"

void init_poly(polynomial_mpz * restrict polynomial)
{
    polynomial->coeffs = calloc(1, sizeof(mpz_t));
    mpz_init_set_ui(polynomial->coeffs[0], 0);

    polynomial->degree = 0;
}

void init_poly_degree(polynomial_mpz * restrict polynomial, const unsigned long degree)
{
    polynomial->coeffs = calloc(degree + 1, sizeof(mpz_t));
    for (size_t i = 0 ; i <= degree ; i++)
    {
        mpz_init_set_ui(polynomial->coeffs[i], 0);
    }

    polynomial->degree = degree;
}

void reduce_polynomial(polynomial_mpz * restrict polynomial) // Delete leading zeros
{
    while (polynomial->degree && !mpz_cmp_ui(polynomial->coeffs[0], 0))
    {
        mpz_t *tmp_array = calloc(polynomial->degree, sizeof(mpz_t));

        for (size_t i = 0 ; i < polynomial->degree ; i++)
        {
            mpz_init_set(tmp_array[i], polynomial->coeffs[i+1]);
        }

        for (size_t i = 0 ; i <= polynomial->degree ; i++)
        {
            mpz_clear(polynomial->coeffs[i]);
        }

        free(polynomial->coeffs);

        polynomial->coeffs = tmp_array;

        polynomial->degree--;
    }
}

void reduce_polynomial_last(polynomial_mpz * restrict polynomial) // Delete last zeros
{
    while (polynomial->degree && !mpz_cmp_ui(polynomial->coeffs[polynomial->degree], 0))
    {
        mpz_t *tmp_array = calloc(polynomial->degree, sizeof(mpz_t));

        for (size_t i = 0 ; i < polynomial->degree ; i++)
        {
            mpz_init_set(tmp_array[i], polynomial->coeffs[i]);
        }

        for (size_t i = 0 ; i <= polynomial->degree ; i++)
        {
            mpz_clear(polynomial->coeffs[i]);
        }

        free(polynomial->coeffs);

        polynomial->coeffs = tmp_array;

        polynomial->degree--;
    }
}

void set_coeff(polynomial_mpz * restrict polynomial, const mpz_t number, const unsigned long index)
{
    if (polynomial->degree >= index)
    {
        mpz_set(polynomial->coeffs[polynomial->degree - index], number);
    }
    else
    {
        mpz_t *new_array = calloc(index + 1, sizeof(mpz_t));

        mpz_init_set(new_array[0], number);

        for (size_t i = 1 ; i < index - polynomial->degree ; i++)
        {
            mpz_init_set_ui(new_array[i], 0);
        }

        for (size_t i = 0 ; i <= polynomial->degree ; i++)
        {
            mpz_init_set(new_array[index - polynomial->degree + i], polynomial->coeffs[i]);
            mpz_clear(polynomial->coeffs[i]);
        }

        free(polynomial->coeffs);

        polynomial->coeffs = new_array;
        polynomial->degree = index;
    }
}

void copy_polynomial(polynomial_mpz *polynomial1, const polynomial_mpz *polynomial2) // Copy polynomial2 into polynomial1
{
    for (size_t i = 0 ; i <= polynomial1->degree ; i++)
    {
        mpz_clear(polynomial1->coeffs[i]);
    }

    free(polynomial1->coeffs);

    polynomial1->coeffs = calloc(polynomial2->degree + 1, sizeof(mpz_t));

    for (size_t i = 0 ; i <= polynomial2->degree ; i++)
    {
        mpz_init_set(polynomial1->coeffs[i], polynomial2->coeffs[i]);
    }

    polynomial1->degree = polynomial2->degree;
}

void reset_polynomial(polynomial_mpz * restrict polynomial)
{
    for (size_t i = 0 ; i <= polynomial->degree ; i++)
    {
        mpz_clear(polynomial->coeffs[i]);
    }

    free(polynomial->coeffs);

    polynomial->coeffs = calloc(1, sizeof(mpz_t));
    mpz_init_set_ui(polynomial->coeffs[0], 0);

    polynomial->degree = 0;
}

void free_polynomial(polynomial_mpz * restrict polynomial)
{
    for (size_t i = 0 ; i <= polynomial->degree ; i++)
    {
        mpz_clear(polynomial->coeffs[i]);
    }

    free(polynomial->coeffs);
    polynomial->coeffs = NULL;

    polynomial->degree = 0;
}

void print_polynomial(const polynomial_mpz * restrict polynomial)
{
    printf("[");

    for (size_t i = 0 ; i < polynomial->degree ; i++)
    {
        gmp_printf("%Zd, ", polynomial->coeffs[i]);
    }

    gmp_printf("%Zd]\n", polynomial->coeffs[polynomial->degree]);
}