#ifndef UTILS_POLYNOMIAL_SELECTION
#define UTILS_POLYNOMIAL_SELECTION

#include <gmp.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "polynomial_functions.h"

void evaluate_legendre(
    const unsigned long n,
    const double x,
    double * restrict P,
    double * restrict dP
);

void nodes_and_weights_gauss(
    const unsigned long n,
    double * restrict x,
    double * restrict w
);

void get_Lnorm(
    mpf_t res,
    const polynomial_mpz * restrict F,
    const unsigned long skew,
    const unsigned long smoothness_bound,
    const mpf_t ln2,
    const mpf_t e
);

double get_Escore(
    const polynomial_mpz * restrict f,
    const polynomial_mpz * restrict g,
    const double alpha,
    const unsigned long smoothness_bound,
    const unsigned long x_limit,
    const unsigned long y_limit,
    dickman_table * restrict table,
    const mpf_t ln2,
    const mpf_t e
);

void minimize_Lnorm(
    polynomial_mpz * restrict f,
    unsigned long skew_factor,
    const unsigned long smoothness_bound,
    mpz_t m0,
    mpz_t m1,
    const mpf_t ln2,
    const mpf_t e
);

void get_sieve_region(
    const polynomial_mpz * restrict f,
    const unsigned long smoothness_bound,
    unsigned long * restrict max_a_norm,
    unsigned long * restrict skew_factor,
    const mpf_t ln2,
    const mpf_t e
);

void build_poly_coeffs(
    polynomial_mpz * restrict f,
    const mpz_t m0,
    const mpz_t m1,
    const mpz_t a_d,
    const mpz_t n,
    const unsigned long d
);

void get_alpha_score(
    double * restrict res,
    polynomial_mpz * restrict f,
    dyn_array_classic * restrict primes,
    const mpf_t ln2,
    const mpf_t e,
    gmp_randstate_t state
);

void compute_m_mu(
    mpz_t res,
    const mpz_t m0,
    size_t l,
    unsigned long d,
    const mpz_t components[l][d],
    const unsigned long * restrict vec
);

void compute_e(
    mpz_t m0,
    unsigned long nb_roots,
    unsigned long d,
    mpz_t roots_used[nb_roots][d],
    mpz_t prod,
    mpz_t a_d,
    mpz_t n,
    mpz_t e[nb_roots][d]
);


#endif // UTILS_POLYNOMIAL_SELECTION