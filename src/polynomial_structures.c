#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "polynomial_structures.h"

void init_poly(polynomial_mpz *polynomial)
{
    polynomial->coeffs = calloc(1, sizeof(mpz_t));
    mpz_init_set_ui(polynomial->coeffs[0], 0);

    polynomial->degree = 0;
}

void init_poly_degree(polynomial_mpz *polynomial, unsigned long degree)
{
    polynomial->coeffs = calloc(degree+1, sizeof(mpz_t));
    for (size_t i = 0 ; i <= degree ; i++)
    {
        mpz_init_set_ui(polynomial->coeffs[i], 0);
    }

    polynomial->degree = degree;
}

void set_coeff(polynomial_mpz *polynomial, mpz_t number, unsigned long index)
{
    if (polynomial->degree >= index)
    {
        mpz_set(polynomial->coeffs[index], number);
    }
    else
    {
        mpz_t *new_array = calloc(index+1, sizeof(mpz_t));

        for (size_t i = 0 ; i <= index ; i++)
        {
            mpz_init(new_array[i]);
        }

        for (size_t i = polynomial->degree+1 ; i <= index ; i++)
        {
            mpz_set_ui(new_array[i], 0);void print_polynomial(polynomial_mpz *polynomial);
        }

        for (size_t i = 0 ; i <= polynomial->degree ; i++)
        {
            mpz_set(new_array[i], polynomial->coeffs[i]);
            mpz_clear(polynomial->coeffs[i]);
        }

        free(polynomial->coeffs);

        polynomial->coeffs = new_array;
        polynomial->degree = index;
        
        mpz_set(polynomial->coeffs[index], number);
    }
}

void copy_polynomial(polynomial_mpz *polynomial1, polynomial_mpz *polynomial2) // Copy polynomial2 into polynomial1
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

void free_polynomial(polynomial_mpz *polynomial)
{
    for (size_t i = 0 ; i <= polynomial->degree ; i++)
    {
        mpz_clear(polynomial->coeffs[i]);
    }

    free(polynomial->coeffs);
    polynomial->coeffs = NULL;

    polynomial->degree = 0;
}

void print_polynomial(polynomial_mpz *polynomial)
{
    printf("[");

    for (size_t i = 0 ; i < polynomial->degree ; i++)
    {
        gmp_printf("%Zd, ", polynomial->coeffs[i]);
    }

    gmp_printf("%Zd]\n", polynomial->coeffs[polynomial->degree]);
}