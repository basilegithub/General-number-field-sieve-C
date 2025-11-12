#include <gmp.h>
#include <stdlib.h>

void compute_e(mpf_t e)
{
    mpf_set_ui(e, 1);
    mpf_set_prec(e, 96);
    mpf_t fac, tmp;
    mpf_init_set_ui(fac, 1);
    mpf_init(tmp);
    for (unsigned long i = 1 ; i < 50 ; i++)
    {
        mpf_mul_ui(fac, fac, i);
        mpf_ui_div(tmp, 1, fac);
        mpf_add(e, e, tmp);
    }
    mpf_clears(fac, tmp, NULL);
}

void initialize_params(gmp_randstate_t state, mpf_t ln10, mpf_t ln2, mpf_t e)
{
    mpz_t n, b;
    mpz_init_set_ui(n, 3);
    mpz_init_set_ui(b, 82939);
    mpz_powm_ui(n, n, rand()*23%103, b);
    gmp_randseed(state, n);

    mpf_init_set_d(ln10,2.302585092994046);
    mpf_init_set_d(ln2,0.6931471805599453);

    compute_e(e);

    mpz_clears(n, b, NULL);
}

void natural_log(mpf_t res, mpf_t x, mpf_t ln2, mpf_t e)
{
    mpz_t tmp;

    mpz_init(tmp);
    mpz_set_f(tmp, x);
    mpz_set_ui(tmp, mpz_sizeinbase(tmp, 2));

    mpf_t a, tmpf;

    mpf_inits(a, tmpf, NULL);
    mpf_set_prec(a, 32);
    mpf_set_z(a, tmp);
    mpf_mul(a, a, ln2);

    for (unsigned long i = 0 ; i < 5 ; i++)
    {
        myexp(tmpf, a, e);
        mpf_div(tmpf, x, tmpf);
        mpf_sub_ui(a, a, 1);
        mpf_add(a, a, tmpf);
    }
    mpf_set(res, a);

    mpz_clear(tmp);
    mpf_clears(tmpf, a, NULL);
}