#include <gmp.h>
#include <stdlib.h>

#include "dynamic_arrays.h"
#include "algebraic_base.h"
#include "NFS_relations.h"
#include "utils.h"

unsigned int compute_degree(mpz_t n, mpf_t ln2, mpf_t e)
{
    mpf_t tmpf, tmpf2, tmpf3;
    mpf_inits(tmpf, tmpf2, tmpf3, NULL);

    mpf_set_z(tmpf3, n);

    natural_log(tmpf, tmpf3, ln2, e);
    natural_log(tmpf2, tmpf, ln2, e);

    mpf_mul_ui(tmpf3, tmpf, 3);

    mpf_div(tmpf, tmpf3, tmpf2);

    nth_root(tmpf2, tmpf, 3);

    unsigned long res = mpf_get_ui(tmpf2);

    mpf_clears(tmpf, tmpf2, tmpf3, NULL);

    return res;
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
    nth_root(tmpf, tmpf, 3);

    natural_log(tmpf2, tmpf2, ln2, e);
    natural_log(tmpf2, tmpf2, ln2, e);
    mpf_pow_ui(tmpf2, tmpf2, 2);
    nth_root(tmpf2, tmpf2, 3);

    mpf_mul(tmpf3, tmpf, tmpf2);

    myexp(tmpf4, tmpf3, e);
    mpf_div_ui(tmpf4, tmpf4, 6); // Small division by constant

    mpz_set_f(smooth_bound, tmpf4);

    mpf_clears(tmpf, tmpf2, tmpf3, tmpf4, NULL);
}

void compute_logs(dyn_array_classic *logs, dyn_array_classic primes)
{
    mpz_t tmp;
    mpz_init(tmp);

    for (size_t i = 0 ; i < primes.len ; i++)
    {
        mpz_set_ui(tmp, primes.start[i]);
        append_classic(logs, mpz_sizeinbase(tmp, 2) - 1);
    }

    mpz_clear(tmp);
}

void compute_free_relations(
    nfs_relations * restrict relations,
    const algebraic_base * restrict alg_base,
    const unsigned long * restrict divide_leading,
    mpz_t leading_coeff,
    mpz_t m1,
    const size_t len_divide_leading,
    const unsigned long degree
)
{
    algebraic_base_prime *alg_prime = alg_base->start;

    mpz_t tmp;
    mpz_init(tmp);

    while (alg_prime != NULL)
    {
        if (alg_prime->roots.len == d)
        {
            init_new_relation(relations, len_divide_leading);

            mpz_set_ui(tmp, alg_prime->prime);
            set_coeff(&relations->rels[relations->len - 1].poly_f, tmp, 0); // f(x) = p

            mpz_mul(tmp, tmp, leading_coeff);
            set_coeff(&relations->rels[relations->len - 1].poly_g, tmp, 0); // g(x) = c_d * p

            if (!(degree&1))
            {
                mpz_set(relations->rels[relations->len - 1].algebraic_norm, tmp); // algebraic_norm = c_d * p
            }
            else
            {
                mpz_set_ui(tmp, alg_prime->prime);
                mpz_pow_ui(tmp, tmp, degree);
                mpz_mul(tmp, tmp, leading_coeff);
                mpz_set(relations->rels[relations->len - 1].algebraic_norm, tmp); // algebraic_norm = c_d * p^d
            }

            mpz_set_ui(tmp, alg_prime->prime);
            mpz_mul(tmp, tmp, m1);
            mpz_set(relations->rels[relations->len - 1].rational_norm, tmp); // rational_norm = m_1 * p
        }

        alg_prime = alg_prime->next;
    }

    mpz_clear(tmp);
}