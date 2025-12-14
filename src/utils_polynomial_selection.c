#include <gmp.h>
#include <math.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "polynomial_functions.h"
#include "utils.h"

void get_Lnorm(
    mpz_t res,
    const polynomial_mpz * restrict F,
    const unsigned long skew,
    const unsigned long smoothness_bound,
    const mpf_t ln2,
    const mpf_t e
)
{
    unsigned long square_root = isqrt(skew);
    unsigned long n = F->degree + 1;

    mpz_t base_X, base_Y;
    mpz_init_set_ui(base_X, smoothness_bound);
    mpz_init_set_ui(base_Y, smoothness_bound);

    mpz_mul_ui(base_X, base_X, square_root);
    mpz_div_2exp(base_X, base_X, 1);

    mpz_div_ui(base_Y, base_Y, square_root);
    mpz_add_ui(base_Y, base_Y, 1);

    mpz_t current_X, current_Y;
    mpz_init_set(current_X, base_X);
    mpz_init_set(current_Y, base_Y);

    mpz_pow_ui(current_X, current_X, n);

    mpz_mul(base_X, base_X, base_X);
    mpz_mul(base_Y, base_Y, base_Y);

    mpz_set_ui(res, 0);

    mpz_t tmp, tmp2;
    mpz_inits(tmp, tmp2, NULL);

    for (size_t i = 0 ; i < n ; i += 2)
    {
        mpz_mul_2exp(tmp, current_X, 1);
        mpz_sub_ui(tmp2, current_X, 1);

        mpz_mul(tmp, tmp, tmp2);
        mpz_mul(tmp, tmp, F->coeffs[i]);

        mpz_set_ui(tmp2, (i+1)*(n-i));
        mpz_div(tmp, tmp2);

        mpz_add(res, res, tmp);

        mpz_divexact(current_X, base_X);
        mpz_mul(current_Y, base_Y);
    }

    mpz_abs(res, res);
    
    mpz_t res_f, x;
    mpf_init(res_f);
    mpf_init(x);

    mpf_set_prec(x, 2048);
    mpf_set_z(x, res);

    natural_log(res_f, x, ln2, e);
    mpz_set_f(res, res_f);

    mpz_clears(base_X, base_Y, current_X, current_Y, tmp, tmp, NULL);
    mpf_clears(res_f, x, NULL);
}