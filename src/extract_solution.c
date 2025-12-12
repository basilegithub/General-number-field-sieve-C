#include <gmp.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "NFS_relations.h"
#include "polynomial_functions.h"
#include "utils.h"
#include "square_root.h"

void extract_solution(
    mpz_t factor,
    const nfs_relations * restrict relations,
    const bool * restrict kernel_vector,
    const dyn_array_classic * restrict rational_primes,
    const polynomial_mpz f_x,
    const polynomial_mpz f_prime_sq,
    const mpz_t leading_coeff,
    const mpz_t n,
    const mpz_t m0,
    const mpz_t m1,
    const mpz_t f_prime_eval,
    const unsigned long inert_prime,
    const unsigned long max_a_size,
    gmp_randstate_t state,
    FILE *logfile
)
{
    mpz_t f_norm;
    mpz_init(f_norm);

    for (size_t i = 0 ; i <= f_x.degree ; i++)
    {
        mpz_addmul(f_norm, f_x.coeffs[i], f_x.coeffs[i]);
    }
    mpz_sqrt(f_norm, f_norm);

    mpz_t fd;
    mpz_init_set_ui(fd, f_x.degree);
    mpz_pow_ui(fd, fd, 3);
    mpz_sqrt(fd, fd);

    mpz_t x, rational_square_root;
    mpz_init_set(x, f_prime_eval);
    mpz_init(rational_square_root);

    extract_rational_square_root(rational_square_root, relations, kernel_vector, n, rational_primes);

    mpz_mul(x, x, rational_square_root);
    mpz_mod(x, x, n);

    unsigned long S = 0;

    polynomial_mpz algebraic_square, tmp_poly, tmp_poly2;
    init_poly(&algebraic_square);
    init_poly(&tmp_poly);
    init_poly(&tmp_poly2);

    copy_polynomial(&algebraic_square, &f_prime_sq);

    for (size_t i = 0 ; i < relations->len ; i++)
    {
        if (kernel_vector[i])
        {
            poly_prod(&tmp_poly, &algebraic_square, &relations->rels[i].poly_g);
            poly_div(&algebraic_square, &tmp_poly, &f_x);

            S += relations->rels[i].nb_relations;
        }
    }

    mpz_t tmp_mpz, tmp_mpz2;
    mpz_inits(tmp_mpz, tmp_mpz2, NULL);

    mpz_powm_ui(tmp_mpz, leading_coeff, S>>1, n);

    mpz_mul(x, x, tmp_mpz);
    mpz_mod(x, x, n);

    // gmp_printf("x = %Zd\n", x);

    mpz_t coeff_bound;
    mpz_init_set_ui(coeff_bound, 1);

    for (size_t i = 0 ; i < f_x.degree ; i++)
    {
        mpz_pow_ui(tmp_mpz, f_norm, f_x.degree - i);
        mpz_mul(tmp_mpz, tmp_mpz, fd);
        mpz_mul_ui(tmp_mpz2, leading_coeff, max_a_size);
        mpz_mul_2exp(tmp_mpz2, tmp_mpz2, 1);
        mpz_mul(tmp_mpz2, tmp_mpz2, f_norm);
        mpz_pow_ui(tmp_mpz2, tmp_mpz2, S>>1);
        mpz_mul(tmp_mpz, tmp_mpz, tmp_mpz2);

        if (mpz_cmp(tmp_mpz, coeff_bound) > 0) mpz_set(coeff_bound, tmp_mpz);
    }

    mpz_t algebraic_square_root;
    mpz_init(algebraic_square_root);

    extract_algebraic_square_root(algebraic_square_root, f_x, algebraic_square, m0, m1, leading_coeff, coeff_bound, inert_prime, state, logfile);
    mpz_mod(algebraic_square_root, algebraic_square_root, n);

    mpz_powm_ui(tmp_mpz, m1, S>>1, n);

    mpz_mul(algebraic_square_root, algebraic_square_root, tmp_mpz);
    mpz_mod(algebraic_square_root, algebraic_square_root, n);

    // gmp_printf("y = %Zd\n", algebraic_square_root);

    mpz_add(factor, x, algebraic_square_root);

    free_polynomial(&tmp_poly);
    free_polynomial(&tmp_poly2);
    free_polynomial(&algebraic_square);

    mpz_clears(f_norm, fd, tmp_mpz, tmp_mpz2, coeff_bound, rational_square_root, x, algebraic_square_root, NULL);
}