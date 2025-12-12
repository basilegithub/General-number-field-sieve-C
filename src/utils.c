#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#include "dynamic_arrays.h"
#include "utils.h"

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

void recursive_exp(mpf_t res, mpz_t pow, const mpf_t e)
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

void myexp(mpf_t res, mpf_t x, const mpf_t e)
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

void natural_log(mpf_t res, mpf_t x, const mpf_t ln2, const mpf_t e)
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

void nth_root(mpf_t r, const mpf_t x, const unsigned long n)
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

    mpf_set(r, tmpf3);

    mpf_clears(tmpf, tmpf2, tmpf3, tmpf4, NULL);
}

unsigned long gcd(unsigned long a, unsigned long b)
{
    while (b != 0) {
        unsigned long t = b;
        b = a % b;
        a = t;
    }
    return a;
}

bool fermat_primality(const mpz_t n)
{
    mpz_t res, base, exponent;
    mpz_init_set_ui(base, 2);
    mpz_inits(res, exponent, NULL);

    mpz_sub_ui(exponent, n, 1);

    mpz_powm(res, base, exponent, n);
    bool to_return = (bool) !mpz_cmp_ui(res, 1);

    mpz_clears(res, base, exponent, NULL);

    return to_return;
}

int my_legendre(const mpz_t n, unsigned long p)
{
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mod_ui(tmp, n, p);

    unsigned long tmpl;
    tmpl = mpz_get_ui(tmp);

    int t = 1;
    unsigned long tmps;
    while (tmpl)
    {
        while (!(tmpl&1))
        {
            tmpl >>= 1;
            if (p%8 == 3 || p%8 == 5) t = -t;
        }
        tmps = tmpl;
        tmpl = p;
        p = tmps;
        if (tmpl%4 == p%4 && p%4 == 3) t = -t;
        tmpl %= p;
    }
    mpz_clear(tmp);
    if (p == 1) return t;
    return 0;
}

void sqrt_mod(mpz_t n, const unsigned long p, gmp_randstate_t state)
{
    mpz_t z, tmp, P_value;
    mpz_inits(z, tmp, NULL);

    unsigned long P;
    unsigned long tmp3;
    unsigned long r = 0;

    mpz_mod_ui(n, n, p);
    P = p-1;
    mpz_init_set_ui(P_value, P);
    mpz_set_ui(tmp, p);

    mpz_urandomm(z, state, tmp);
    if (mpz_cmp_ui(z, 1) != 1) mpz_set_ui(z, 2);

    while (my_legendre(z, p) != -1)
    {
        mpz_urandomm(z, state, tmp);
        if (mpz_cmp_ui(z, 1) != 1) mpz_set_ui(z, 2);
    }

    while (!(P&1))
    {
        P >>= 1;
        r++;
    }
    mpz_t generator, lambda, omega, res, m, two_mpz;
    mpz_inits(generator, lambda, omega, NULL);

    mpz_powm_ui(generator, z, P, tmp);
    mpz_powm_ui(lambda, n, P, tmp);

    tmp3 = (P+1)>>1;
    mpz_powm_ui(omega, n, tmp3, tmp);

    mpz_init_set_ui(res, 0);
    mpz_init(m);
    mpz_init_set_ui(two_mpz, 2);

    mpz_t tmp_l, tmp2;
    mpz_inits(tmp_l, tmp2, NULL);

    while (1)
    {
        if (!mpz_cmp_ui(lambda, 0))
        {
            mpz_set_ui(n, 0);
            break;
        }

        if (!mpz_cmp_ui(lambda, 1))
        {
            mpz_set(n, omega);
            break;
        }

        mpz_set_ui(m, 1);
        while (mpz_cmp_ui(m, r) < 0)
        {
            mpz_powm(tmp_l, two_mpz, m, P_value);
            mpz_powm(tmp_l, lambda, tmp_l, tmp);

            if (!mpz_cmp_ui(tmp_l, 1)) break;

            mpz_add_ui(m, m, 1);
        }

        mpz_ui_sub(tmp2, r-1, m);
        mpz_powm(tmp_l, two_mpz, tmp2, P_value);
        mpz_mul_2exp(tmp2, tmp_l, 1);

        mpz_powm(tmp2, generator, tmp2, tmp);
        mpz_mul(lambda, lambda, tmp2);
        mpz_mod(lambda, lambda, tmp);

        mpz_powm(tmp2, generator, tmp_l, tmp);
        mpz_mul(omega, omega, tmp2);
        mpz_mod(omega, omega, tmp);
    }

    mpz_clears(z, tmp, tmp2, P_value, generator, lambda, omega, res, m, two_mpz, NULL);
}

void convert_to_vec(mpz_t embedding, const unsigned long relations_len, bool * restrict tmp_vec)
{
    mpz_t tmp, tmp2;
    mpz_inits(tmp, tmp2, NULL);
    mpz_set_ui(tmp2, 1);

    for (size_t i = 0 ; i < relations_len ; i++)
    {
        mpz_and(tmp, embedding, tmp2);
        tmp_vec[relations_len-i-1] = (bool) mpz_get_ui(tmp);

        mpz_div_2exp(embedding, embedding, 1);
    }

    mpz_clears(tmp, tmp2, NULL);
}