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
    const size_t l,
    const unsigned long d,
    const mpz_t components[l][d],
    const unsigned long * restrict vec
);

void compute_e_array(
    const mpz_t m0,
    const unsigned long nb_roots,
    const unsigned long d,
    const mpz_t roots_used[nb_roots][d],
    const mpz_t prod,
    const mpz_t a_d,
    const mpz_t n,
    mpz_t e[nb_roots][d]
);

void compute_f(
    const mpz_t n,
    const mpz_t a_d,
    const mpz_t m0,
    const unsigned long d,
    const mpz_t prod,
    const unsigned long nb_roots,
    const mpz_t roots_used[nb_roots][d],
    const mpz_t e[nb_roots][d],
    mpf_t f[nb_roots][d],
    mpf_t f0
);

void create_first_array(
    const unsigned long nb_roots,
    const unsigned long d,
    const mpf_t f0,
    const mpf_t f[nb_roots][d],
    const unsigned long nb_rows,
    const unsigned long nb_cols,
    mpf_t * restrict array1_u_values,
    unsigned long array1_indices[nb_rows][nb_cols],
);

void create_second_array(
    const unsigned long nb_roots,
    const unsigned long vec_len,
    const unsigned long d,
    const unsigned long nb_rows,
    const unsigned long nb_cols,
    const mpf_t f[nb_roots][d],
    mpf_t * restrict array2_u_values,
    unsigned long array2_indices[nb_rows][nb_cols]
);


#endif // UTILS_POLYNOMIAL_SELECTION