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

    // Initialize random seed

    srand(time(NULL));
    mpz_t random_a, random_b;
    mpz_init_set_ui(random_a, 3);
    mpz_init_set_ui(random_b, 82939);
    mpz_powm_ui(random_a, random_a, rand()*23%103, random_b);

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed(state, random_a);

    // Interactive input

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