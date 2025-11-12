#include <gmp.h>
#include <stdlib.h>

#include "utils.h"

unsigned int compute_degree(mpz_t n, mpf_t ln2, mpf_t e)
{
    mpf_t tmp, tmp2, tmp3;
    mpf_inits(tmp, tmp2, tmp3, NULL);

    mpf_set_z(tmp3, n);

    natural_log(tmp, tmp3, ln2, e);
    natural_log(tmp2, tmp, ln2, e);

    mpf_mul_ui(tmp3, tmp, 3);

    mpf_div(tmp, tmp3, tmp2);

    unsigned long res = (mpf_get_ui(tmp))^(1/3);

    mpz_clears(tmp, tmp2, tmp3, NULL);

    return res;
}