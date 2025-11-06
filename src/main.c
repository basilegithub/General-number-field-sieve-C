// Main file of the project

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <gmp.h>

int main()
{
    printf("program started\n");

    char number[100];
    printf("enter the number you wish to factor : ");
    fgets(number, 100, stdin);

    mpz_t n;

    mpz_init(n);

    mpz_init_set_str(n, number, 10);

    // The GNFS algorithm is split in 4 main steps

    // Polynomial selection

    // Sieving

    // Linear algebra

    // Square root extraction



    mpz_clear(n);
}