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

void recursive_exp(mpf_t res, mpz_t pow, mpf_t e)
{
    if (!mpz_cmp_ui(pow, 0)) mpf_set_ui(res, 1);
    else if (!mpz_cmp_ui(pow, 1)) mpf_set(res, e);
    else
    {
        mpf_t tmpf;
        mpf_init(tmpf);
        mpz_t tmp;
        mpz_init(tmp);

        if (mpz_even_p(pow))
        {
            mpz_div_ui(tmp, pow, 2);
            recursive_exp(tmpf, tmp, e);
            mpf_mul(res, tmpf, tmpf);
        }
        else
        {
            mpz_div_ui(tmp, pow, 2);
            recursive_exp(tmpf, tmp, e);
            mpf_mul(tmpf, tmpf, tmpf);
            mpf_mul(res, e, tmpf);
        }
        mpf_clear(tmpf);
        mpz_clear(tmp);
    }
}

void myexp(mpf_t res, mpf_t x, mpf_t e)
{
    int flag = 0;
    if (mpf_cmp_ui(x, 0) < 0)
    {
        flag = 1;
        mpf_neg(x, x);
    }
    mpf_t tmp_res, rest;
    mpz_t tmp;

    mpf_inits(tmp_res, rest, NULL);
    mpz_init(tmp);

    mpz_set_f(tmp, x);
    mpf_set_z(rest, tmp);
    mpf_sub(rest, x, rest);

    recursive_exp(tmp_res, tmp, e);

    mpf_t tmp_rest, fact, tmpf, tmpf2;

    mpf_init_set(tmp_rest, rest);
    mpf_init_set_ui(fact, 1);
    mpf_init_set_ui(tmpf, 1);
    mpf_init(tmpf2);

    for (unsigned long i = 2 ; i <= 10 ; i++)
    {
        mpf_div(tmpf2, tmp_rest, fact);
        mpf_add(tmpf, tmpf, tmpf2);
        mpf_mul_ui(fact, fact, i);
        mpf_mul(tmp_rest, tmp_rest, rest);
    }
    mpf_mul(res, tmp_res, tmpf);
    if (flag)
    {
        mpf_set_ui(tmpf, 1);
        mpf_div(res, tmpf, res);
    }
    mpz_clear(tmp);
    mpf_clears(tmpf, tmpf2, rest, tmp_res, tmp_rest, fact, NULL);
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

    for (size_t i = 0 ; i < 100 ; i++)
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