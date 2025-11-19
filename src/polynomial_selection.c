#include <gmp.h>

#include "polynomial_structures.h"
#include "utils.h"

void basic_polynomial_selection(polynomial_mpz *polynomial, mpz_t n, mpz_t m0, mpz_t m1, unsigned long d)
{
    mpz_t d_root, tmp, tmp2, tmp3;
    mpz_inits(d_root, tmp, tmp2, tmp3, NULL);

    mpf_t tmpf;
    mpf_init(tmpf);
    mpf_set_z(tmpf, n);
    nth_root(tmpf, tmpf, d);

    mpz_set_f(d_root, tmpf);

    mpz_set(m0, d_root);
    mpz_set_ui(m1, 1);

    mpf_clear(tmpf);

    mpz_set(tmp, n);

    for (unsigned long i = 0 ; i <= d ; i++)
    {
        mpz_pow_ui(tmp2, d_root, d-i);

        mpz_div(tmp3, tmp, tmp2);

        set_coeff(polynomial, tmp3, d-i);

        mpz_mod(tmp, tmp, tmp2);
    }

    mpz_clears(d_root, tmp, tmp2, tmp3, NULL);
}