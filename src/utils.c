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

void nth_root(mpf_t r, mpf_t x, unsigned long n)
{
    mpf_t tmpf, tmpf2, tmpf3, tmpf4;
    mpf_inits(tmpf, tmpf2, tmpf3, tmpf4, NULL);

    mpf_set_ui(tmpf, 1);
    mpf_set(tmpf2, x);

    mpf_add(tmpf3, tmpf, tmpf2);
    mpf_div_2exp(tmpf3, tmpf3, 1);

    mpf_pow_ui(tmpf4, tmpf3, n);

    for (size_t i = 0 ; i < 50 ; i++)
    {
        if (mpf_cmp(tmpf4, x) < 0) mpf_set(tmpf, tmpf3);
        else mpf_set(tmpf2, tmpf3);

        mpf_add(tmpf3, tmpf, tmpf2);
        mpf_div_2exp(tmpf3, tmpf3, 1);

        mpf_pow_ui(tmpf4, tmpf3, n);
    }

    mpf_set(r, tmpf4);

    mpf_clears(tmpf, tmpf2, tmpf3, tmpf4, NULL);
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


void compute_smooth_bound(mpz_t n, mpz_t smooth_bound, mpf_t ln2, mpf_t e)
{
    mpf_t tmpf, tmpf2, tmpf3, tmpf4;
    mpf_inits(tmpf, tmpf2, tmpf3, tmpf4, NULL);

    mpf_set_z(tmpf, n);
    mpf_set_z(tmpf2, n);

    natural_log(tmpf, tmpf, ln2, e);
    mpf_mul_ui(tmpf, tmpf, 8);
    mpf_div_ui(tmpf, tmpf, 9);
    

    natural_log(tmpf2, tmpf2, ln2, e);
    natural_log(tmpf2, tmpf2, ln2, e);

    mpf_clears(tmpf, tmpf2, tmpf3, tmpf4, NULL);
}