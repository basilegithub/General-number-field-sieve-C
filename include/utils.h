#ifndef UTILS_H
#define UTILS_H

#include <gmp.h>
#include <stdbool.h>

#include "dynamic_arrays.h"

void compute_e(mpf_t e);
void recursive_exp(mpf_t res, mpz_t pow, mpf_t e);
void myexp(mpf_t res, mpf_t x, mpf_t e);
void natural_log(mpf_t res, mpf_t x, mpf_t ln2, mpf_t e);
void nth_root(mpf_t r, mpf_t x, unsigned long n);
void initialize_params(gmp_randstate_t state, mpf_t ln10, mpf_t ln2, mpf_t e);
void compute_smooth_bound(mpz_t n, mpz_t smooth_bound, mpf_t ln2, mpf_t e);
unsigned long gcd(unsigned long a, unsigned long b);
bool fermat_primality(mpz_t n);
int my_legendre(mpz_t n, unsigned long p);
void sqrt_mod(mpz_t n, const unsigned long p, gmp_randstate_t state);
void convert_to_vec(mpz_t embedding, unsigned long relations_len, bool tmp_vec[relations_len]);

#endif // UTILS_H