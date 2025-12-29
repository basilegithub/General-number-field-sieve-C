#ifndef UTILS_H
#define UTILS_H

#include <gmp.h>
#include <stdbool.h>

#include "dynamic_arrays.h"

void compute_e(mpf_t e);

void recursive_exp(mpf_t res, mpz_t pow, const mpf_t e);

void myexp(mpf_t res, mpf_t x, const mpf_t e);

void natural_log(mpf_t res, mpf_t x, const mpf_t ln2, const mpf_t e);

unsigned long isqrt(unsigned long n);

void nth_root(mpf_t r, const mpf_t x, const unsigned long n);

unsigned long my_power(unsigned long a, unsigned long n);

unsigned long gcd(unsigned long a, unsigned long b);

bool fermat_primality(const mpz_t n);

int my_legendre(const mpz_t n, unsigned long p);

void sqrt_mod(mpz_t n, const unsigned long p, gmp_randstate_t state);

void my_factorial(mpz_t res, unsigned long n);

void gamma_function_div2(
    mpf_t res,
    unsigned long k
);

void central(
    mpf_t res,
    unsigned long k,
    double x,
    const mpf_t e
);

void non_central(
    mpf_t res,
    unsigned long k,
    double l,
    double x,
    const mpf_t e
);

void binom_coeff(
    mpz_t res,
    const unsigned long k,
    const unsigned long n
);

void convert_to_vec(mpz_t embedding, unsigned long relations_len, bool * restrict tmp_vec);

#endif // UTILS_H