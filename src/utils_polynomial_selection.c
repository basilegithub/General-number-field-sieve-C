#include <gmp.h>
#include <math.h>
#include <stdlib.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "polynomial_functions.h"
#include "dickman_table.h"
#include "utils.h"

void evaluate_legendre(
    const unsigned long n,
    const double x,
    double * restrict P,
    double * restrict dP
)
{
    double P0 = 1.0;
    double P1 = x;
    double Pn;

    for (int k = 2; k <= n; k++) {
        Pn = ((2.0*k - 1.0)*x*P1 - (k - 1.0)*P0) / k;
        P0 = P1;
        P1 = Pn;
    }

    *P = (n == 0) ? P0 : (n == 1) ? P1 : Pn;
    *dP = n * (x * (*P) - P0) / (x*x - 1.0);
}

void nodes_and_weights_gauss(
    const unsigned long n,
    double * restrict x,
    double * restrict w
)
{
    const double EPS = 1e-14;
    int m = (n + 1) / 2;

    for (int i = 0; i < m; i++) {
        // Initial guess (good approximation)
        double z = cos(M_PI * (i + 0.75) / (n + 0.5));
        double z1;

        // Newton-Raphson
        do {
            double P, dP;
            evaluate_legendre(n, z, &P, &dP);
            z1 = z;
            z = z1 - P / dP;
        } while (fabs(z - z1) > EPS);

        // Store symmetric roots
        x[i] = -z;
        x[n - 1 - i] = z;

        // Weight
        double P, dP;
        evaluate_legendre(n, z, &P, &dP);
        double wi = 2.0 / ((1.0 - z*z) * dP * dP);
        w[i] = wi;
        w[n - 1 - i] = wi;
    }
}

void get_Lnorm(
    mpf_t res,
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

    mpz_t tmp, tmp2, res_z;
    mpz_inits(tmp, tmp2, res_z, NULL);

    mpf_t tmp_mpf, tmp_mpf2;
    mpf_inits(tmp_mpf, tmp_mpf2, NULL);

    mpz_set_ui(res_z, 0);

    for (size_t i = 0 ; i < n ; i += 2)
    {
        mpz_mul_2exp(tmp, current_X, 1);
        mpz_sub_ui(tmp2, current_X, 1);

        mpz_mul(tmp, tmp, tmp2);
        mpz_mul(tmp, tmp, F->coeffs[i]);

        mpz_set_ui(tmp2, (i+1)*(n-i));
        mpz_div(tmp, tmp, tmp2);

        mpz_add(res_z, res_z, tmp);

        mpz_divexact(current_X, current_X, base_X);
        mpz_mul(current_Y, current_Y, base_Y);
    }

    mpz_abs(res_z, res_z);
    
    mpf_t x;
    mpf_init(x);

    mpf_set_prec(x, 2048);
    mpf_set_z(x, res_z);

    natural_log(res, x, ln2, e);

    mpz_clears(base_X, base_Y, current_X, current_Y, tmp, tmp, NULL);
    mpf_clears(x, tmp_mpf, tmp_mpf2, NULL);
}

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
)
{
    unsigned long n = 32;

    double log_B;

    mpf_t tmpf, tmpf2;
    mpf_init_set_ui(tmpf, smoothness_bound);
    mpf_init(tmpf2);

    natural_log(tmpf2, tmpf, ln2, e);

    log_B = mpf_get_d(tmpf2);


    double *weights = calloc(n, sizeof(double));
    double *nodes = calloc(n, sizeof(double));
    nodes_and_weights_gauss(n, nodes, weights);

    double a, b;
    a = 0.0;
    b = 3.141592; // pi

    for (size_t i = 0 ; i < n ; i++)
    {
        // Shift nodes and weights from [-1, 1] to the required interval [0, pi]
        nodes[i] = 0.5 * (b - a) * nodes[i] + 0.5 * (b + a);
        weights[i] = 0.5 * (b - a) * weights[i];
    }

    double res = 0;

    mpz_t X_mpz, Y_mpz, res1, res2;
    mpz_inits(X_mpz, Y_mpz, res1, res2, NULL);

    for (size_t i = 0 ; i < n ; i++)
    {
        double angle = nodes[i];
        double X = (double)x_limit * cos(angle);
        double Y = (double)(y_limit - 1) * sin(angle) + 1;

        mpz_set_d(X_mpz, X);
        mpz_set_d(Y_mpz, Y);

        evaluate_homogeneous(res1, f, Y_mpz, Y_mpz);
        mpz_abs(res1, res1);
        mpf_set_z(tmpf, res1);
        natural_log(tmpf2, tmpf, ln2, e);
        mpf_set_d(tmpf, alpha);
        mpf_sub(tmpf2, tmpf2, tmpf);
        mpf_set_d(tmpf, log_B);
        mpf_div(tmpf2, tmpf2, tmpf);

        mpz_set_d(res1, evaluate_dickman(table, mpf_get_d(tmpf2)));

        evaluate_homogeneous(res2, g, Y_mpz, Y_mpz);
        mpz_abs(res1, res1);
        mpf_set_z(tmpf, res1);
        natural_log(tmpf2, tmpf, ln2, e);
        mpf_set_d(tmpf, log_B);
        mpf_div(tmpf2, tmpf2, tmpf);

        mpz_set_d(res2, evaluate_dickman(table, mpf_get_d(tmpf2)));

        mpz_mul(res1, res1, res2);
        mpf_set_z(tmpf, res1);
        mpf_set_d(tmpf2, weights[i]);
        mpf_mul(tmpf, tmpf, tmpf2);
        
        res += mpf_get_d(tmpf);
    }

    mpz_clears(X_mpz, Y_mpz, res1, res2, NULL);

    mpf_clears(tmpf, tmpf2, NULL);

}

void get_Epscore(
    const polynomial_mpz * restrict f,
    const polynomial_mpz * restrict g,
    const double *alpha,
    const unsigned long smoothness_bound,
    const unsigned long x_limit,
    const unsigned long y_limit,
    dickman_table * restrict table,
    const mpf_t ln2,
    const mpf_t e
)
{
    double k = *alpha;
    double l = alpha[1];
    double c = alpha[2];

    unsigned long n = 32;

    double log_B;

    mpf_t tmpf, tmpf2;
    mpf_init_set_ui(tmpf, smoothness_bound);
    mpf_init(tmpf2);

    natural_log(tmpf2, tmpf, ln2, e);

    log_B = mpf_get_d(tmpf2);

    double *weights_mu = calloc(n, sizeof(double));
    double *nodes_mu = calloc(n, sizeof(double));
    nodes_and_weights_gauss(n, nodes_mu, weights_mu);

    double a_mu, b_mu;
    a_mu = c;
    b_mu = 6*c; // pi

    for (size_t i = 0 ; i < n ; i++)
    {
        // Shift nodes and weights from [-1, 1] to the required interval [c, 6c]
        nodes_mu[i] = 0.5 * (b_mu - a_mu) * nodes_mu[i] + 0.5 * (b_mu + a_mu);
        weights_mu[i] = 0.5 * (b_mu - a_mu) * weights_mu[i];
    }

    double *weights_th = calloc(n, sizeof(double));
    double *nodes_th = calloc(n, sizeof(double));
    nodes_and_weights_gauss(n, nodes_th, weights_th);

    double a_th, b_th;
    a_th = 0;
    b_th = 3.141592; // pi

    for (size_t i = 0 ; i < n ; i++)
    {
        // Shift nodes and weights from [-1, 1] to the required interval [0, pi]
        nodes_th[i] = 0.5 * (b_th - a_th) * nodes_th[i] + 0.5 * (b_th + a_th);
        weights_th[i] = 0.5 * (b_th - a_th) * weights_th[i];
    }

    double res = 0.0;

    mpz_t X_mpz, Y_mpz, res1, res2;
    mpz_inits(X_mpz, Y_mpz, res1, res2, NULL);

    mpf_t non_central_eval;
    mpf_init(non_central_eval);

    for (size_t i = 0 ; i < n ; i++)
    {
        non_central(non_central_eval, k, l, nodes_mu[i], e);

        for (size_t j = 0 ; j < n ; j++)
        {
            double angle = nodes_th[j];

            double X = (double)x_limit * cos(angle);
            double Y = (double)(y_limit - 1) * sin(angle) + 1;

            mpz_set_d(X_mpz, X);
            mpz_set_d(Y_mpz, Y);

            evaluate_homogeneous(res1, f, Y_mpz, Y_mpz);
            mpz_abs(res1, res1);
            mpf_set_z(tmpf, res1);
            natural_log(tmpf2, tmpf, ln2, e);
            mpf_set_d(tmpf, nodes_mu[i] - c);
            mpf_sub(tmpf2, tmpf2, tmpf);
            mpf_set_d(tmpf, log_B);
            mpf_div(tmpf2, tmpf2, tmpf);

            mpz_set_d(res1, evaluate_dickman(table, mpf_get_d(tmpf2)));

            evaluate_homogeneous(res2, g, Y_mpz, Y_mpz);
            mpz_abs(res1, res1);
            mpf_set_z(tmpf, res1);
            natural_log(tmpf2, tmpf, ln2, e);
            mpf_set_d(tmpf, log_B);
            mpf_div(tmpf2, tmpf2, tmpf);

            mpz_set_d(res2, evaluate_dickman(table, mpf_get_d(tmpf2)));

            mpz_mul(res1, res1, res2);
            mpf_set_z(tmpf, res1);
            mpf_set_d(tmpf2, weights_mu[i]);
            mpf_mul(tmpf, tmpf, tmpf2);
            mpf_set_d(tmpf2, weights_th[j]);
            mpf_mul(tmpf, tmpf, tmpf2);
            mpf_set(tmpf2, non_central_eval);
            mpf_mul(tmpf, tmpf, tmpf2);
            
            res += mpf_get_d(tmpf);
        }
    }

    mpz_clears(X_mpz, Y_mpz, res1, res2, NULL);

    mpf_clears(tmpf, tmpf2, non_central_eval, NULL);

}

void minimize_Lnorm(
    polynomial_mpz * restrict f,
    unsigned long skew_factor,
    const unsigned long smoothness_bound,
    mpz_t m0,
    mpz_t m1,
    const mpf_t ln2,
    const mpf_t e
)
{
    polynomial_mpz tmp_poly, tmp_poly2;
    init_poly(&tmp_poly);
    init_poly(&tmp_poly2);

    mpf_t res, tmp_mpf;
    mpf_inits(res, tmp_mpf, NULL);


    size_t k = 1;

    poly_prod(&tmp_poly, f, f);
    get_Lnorm(res, &tmp_poly, skew_factor, smoothness_bound, ln2, e);

    while (k > 0)
    {
        bool flag = false;

        shift(&tmp_poly2, f, -k);
        poly_prod(&tmp_poly, &tmp_poly2, &tmp_poly2);
        get_Lnorm(tmp_mpf, &tmp_poly, skew_factor, smoothness_bound, ln2, e);

        if (mpf_cmp(tmp_mpf, res) < 0)
        {
            mpf_set(res, tmp_mpf);
            copy_polynomial(f, &tmp_poly2);
            flag = true;
            mpz_addmul_ui(m0, m1, k);
            k <<= 1;
        }

        shift(&tmp_poly2, f, k);
        poly_prod(&tmp_poly, &tmp_poly2, &tmp_poly2);
        get_Lnorm(tmp_mpf, &tmp_poly, skew_factor, smoothness_bound, ln2, e);

        if (mpf_cmp(tmp_mpf, res) < 0)
        {
            mpf_set(res, tmp_mpf);
            copy_polynomial(f, &tmp_poly2);
            flag = true;
            mpz_submul_ui(m0, m1, k);
            k <<= 1;
        }

        if (!flag) k >>= 1;
    }

    free_polynomial(&tmp_poly);
    free_polynomial(&tmp_poly2);

    mpf_clears(res, tmp_mpf, NULL);
}

void get_sieve_region(
    const polynomial_mpz * restrict f,
    const unsigned long smoothness_bound,
    unsigned long * restrict max_a_norm,
    unsigned long * restrict skew_factor,
    const mpf_t ln2,
    const mpf_t e
)
{
    polynomial_mpz F;
    init_poly(&F);

    poly_prod(&F, f, f);

    mpz_t tmp_mpz, tmp_mpz2;
    mpz_inits(tmp_mpz, tmp_mpz2, NULL);

    mpf_t ratio_mean, tmp_mpf, tmp_mpf2;
    mpf_inits(ratio_mean, tmp_mpf, tmp_mpf2, NULL);

    mpf_set_ui(ratio_mean, 0);

    for (size_t i = 0 ; i < f->degree ; i++)
    {
        mpf_set_z(tmp_mpf, f->coeffs[i]);
        mpf_abs(tmp_mpf, tmp_mpf);
        mpf_add_ui(tmp_mpf, tmp_mpf, 1);

        mpf_set_z(tmp_mpf2, f->coeffs[i+1]);
        mpf_abs(tmp_mpf2, tmp_mpf2);

        mpf_div(tmp_mpf, tmp_mpf2, tmp_mpf);
        
        mpf_set_d(tmp_mpf2, 1e-7);
        mpf_add(tmp_mpf, tmp_mpf, tmp_mpf2);

        natural_log(tmp_mpf2, tmp_mpf, ln2, e);

        mpf_add(ratio_mean, ratio_mean, tmp_mpf2);
    }

    mpf_div_ui(tmp_mpf, ratio_mean, f->degree);
    myexp(ratio_mean, tmp_mpf, e);

    *skew_factor = 2 * mpf_get_ui(ratio_mean);

    unsigned long k = 1;

    mpf_t best_norm;
    mpf_init_set_si(best_norm, -1);

    while (k > 0)
    {
        bool updated = false;

        if (*skew_factor > k)
        {
            get_Lnorm(tmp_mpf, &F, *skew_factor - k, smoothness_bound, ln2, e);

            if (mpf_cmp_si(best_norm, -1) || mpf_cmp(tmp_mpf, best_norm) < 0)
            {
                *skew_factor = *skew_factor - k;
                mpf_set(best_norm, tmp_mpf);
                updated = true;
                k <<= 1;
            }
        }

        if (*skew_factor + k < smoothness_bound)
        {
            get_Lnorm(tmp_mpf, &F, *skew_factor + k, smoothness_bound, ln2, e);

            if (mpf_cmp_si(best_norm, -1) || mpf_cmp(tmp_mpf, best_norm) < 0)
            {
                *skew_factor = *skew_factor + k;
                mpf_set(best_norm, tmp_mpf);
                updated = true;
                k <<= 1;
            }
        }

        if (!updated) k >>= 1;
    }

    mpf_set_ui(tmp_mpf2, *skew_factor);
    nth_root(tmp_mpf, tmp_mpf2, 2);
    mpf_mul_ui(tmp_mpf, tmp_mpf, smoothness_bound);

    *max_a_norm = mpf_get_ui(tmp_mpf);

    mpz_clears(tmp_mpz, tmp_mpz2, NULL);
    mpf_clears(ratio_mean, tmp_mpf, tmp_mpf2, best_norm, NULL);
    free_polynomial(&F);
}

void build_poly_coeffs(
    polynomial_mpz * restrict f,
    const mpz_t m0,
    const mpz_t m1,
    const mpz_t a_d,
    const mpz_t n,
    const unsigned long d
)
{
    set_coeff(f, a_d, d);
    f->degree = d;

    mpz_t r, delta, tmp_mpz, tmp_mpz2;
    mpz_inits(delta, tmp_mpz, tmp_mpz2, NULL);
    mpz_init_set(r, n);

    for (size_t i = 0 ; i < d ; i++)
    {
        mpz_pow_ui(tmp_mpz, m0, d - i);
        mpz_mul(tmp_mpz, tmp_mpz, f->coeffs[i]);
        mpz_divexact(tmp_mpz, tmp_mpz, m1);

        mpz_sub(r, r, tmp_mpz);

        mpz_set(tmp_mpz, m1);
        mpz_pow_ui(tmp_mpz2, m0, d - 1 - i);
        mpz_invert(tmp_mpz, tmp_mpz, tmp_mpz2);
        mpz_mul(tmp_mpz, tmp_mpz, m1);
        mpz_mul(tmp_mpz, tmp_mpz, r);
        mpz_neg(tmp_mpz, tmp_mpz);

        mpz_mul(tmp_mpz2, tmp_mpz2, m1);
        mpz_mod(delta, tmp_mpz, tmp_mpz2);

        mpz_add(tmp_mpz, r, delta);
        mpz_pow_ui(tmp_mpz2, m0, d - 1 - i);

        mpz_divexact(tmp_mpz, tmp_mpz, tmp_mpz2);
        set_coeff(f, tmp_mpz, d - 1 - i);
    }

    mpz_div_2exp(tmp_mpz, m0, 1);
    mpz_neg(tmp_mpz2, tmp_mpz);

    for (size_t i = 1 ; i <= f->degree ; i++)
    {
        if (mpz_cmp(f->coeffs[i], tmp_mpz) > 0)
        {
            mpz_sub(f->coeffs[i], f->coeffs[i], m0);
            mpz_add(f->coeffs[i - 1], f->coeffs[i - 1], m1);
        }
        else if (mpz_cmp(f->coeffs[i], tmp_mpz2) < 0)
        {
            mpz_add(f->coeffs[i], f->coeffs[i], m0);
            mpz_sub(f->coeffs[i - 1], f->coeffs[i - 1], m1);
        }
    }

    mpz_clears(r, delta, tmp_mpz, tmp_mpz2, NULL);
}

void get_alpha_score(
    double * restrict res,
    polynomial_mpz * restrict f,
    dyn_array_classic * restrict primes,
    const mpf_t ln2,
    const mpf_t e,
    gmp_randstate_t state
)
{
    double E = 0.0;
    double F = 0.0;

    polynomial_mpz f_derivative;
    init_poly(&f_derivative);

    poly_derivative(&f_derivative, f);

    double baseline_term = 0.0;

    mpf_t tmp_mpf, tmp_mpf2;
    mpf_inits(tmp_mpf, tmp_mpf2, NULL);

    mpz_t tmp_mpz, tmp_mpz2, tmp_mpz3;
    mpz_inits(tmp_mpz, tmp_mpz2, tmp_mpz3, NULL);

    size_t upto = f->degree + 5;

    for (size_t i = 0 ; i < primes->len ; i++)
    {
        unsigned long p = primes->start[i];

        mpf_set_ui(tmp_mpf, p);
        natural_log(tmp_mpf2, tmp_mpf, ln2, e);

        double log_p = mpf_get_d(tmp_mpf2);

        baseline_term += log_p/(p+1);

        double Ep = 0.0;
        double Fp = 0.0;

        dyn_array_classic roots;
        init_classic(&roots);

        dyn_array ramified_roots;
        init(&ramified_roots);

        find_roots(f, &roots, p, state);

        for (size_t j = 0 ; j < roots.len ; j++)
        {
            if (!evaluate_mod_p(&f_derivative, roots.start[j], p))
            {
                mpz_set_ui(tmp_mpz, roots.start[j]);
                append(&ramified_roots, tmp_mpz);
                Ep += 1/(double)p;
                Fp += 1/(double)p;
            }
            else
            {
                Ep += 1/(double)(p - 1);
                Fp += (p+1)/(double)((p - 1)*(p-1));
            }
        }

        if (ramified_roots.len)
        {
            mpz_set_ui(tmp_mpz, p);

            for (size_t j = 2 ; j < upto ; j++)
            {
                dyn_array new;
                init(&new);

                for (size_t r = 0 ; r < ramified_roots.len ; r++)
                {
                    mpz_mul_ui(tmp_mpz2, tmp_mpz, p);

                    evaluate_mod_p_mpz(f, tmp_mpz3, ramified_roots.start[r], tmp_mpz2);

                    if (!mpz_cmp_ui(tmp_mpz3, 0))
                    {
                        mpz_set(tmp_mpz2, ramified_roots.start[r]);
                        for (size_t k = 0 ; k < p ; k++)
                        {
                            append(&new, tmp_mpz2);
                            mpz_add(tmp_mpz2, tmp_mpz2, tmp_mpz);
                        }

                        mpf_set_z(tmp_mpf, tmp_mpz);
                        mpf_set_d(tmp_mpf2, 1.0);
                        mpf_div(tmp_mpf, tmp_mpf2, tmp_mpf);

                        Ep += mpf_get_d(tmp_mpf);
                        Fp += (2*j - 1)*mpf_get_d(tmp_mpf);
                    }
                }

                mpz_mul_ui(tmp_mpz, tmp_mpz, p);
                copy_dyn(&ramified_roots, &new);

                free_dyn_array(&new);
            }

            mpz_set_ui(tmp_mpz2, p);
            mpz_pow_ui(tmp_mpz2, tmp_mpz2, upto - 2);
            mpz_mul_ui(tmp_mpz2, tmp_mpz2, p - 1);

            mpf_set_z(tmp_mpf, tmp_mpz2);
            mpf_set_ui(tmp_mpf2, ramified_roots.len);
            mpf_div(tmp_mpf, tmp_mpf2, tmp_mpf);
            Ep += mpf_get_d(tmp_mpf);

            mpf_set_d(tmp_mpf, (double)upto*upto + 2*upto/(double)(p-1) + (p+1)/(double)((p-1)*(p-1)));
            mpf_set_z(tmp_mpf2, tmp_mpz);
            mpf_div(tmp_mpf, tmp_mpf, tmp_mpf2);
            mpf_mul_ui(tmp_mpf, tmp_mpf, ramified_roots.len);

            Fp += mpf_get_d(tmp_mpf);
        }

        if (mpz_divisible_ui_p(f->coeffs[0], p))
        {
            polynomial_mpz frev;
            init_poly(&frev);
            for (size_t i = 0 ; i <= f->degree ; i++)
            {
                set_coeff(&frev, f->coeffs[i], i);
            }

            polynomial_mpz frev_derivative;
            init_poly(&frev_derivative);
            poly_derivative(&frev_derivative, &frev);

            if (evaluate_mod_p(&frev_derivative, 0, p))
            {
                Ep += 1/(double)(p - 1);
                Fp += (p + 1)/(double)((p - 1) * (p - 1));
            }
            else
            {
                Ep += 1/(double)p;
                Fp += 1/(double)p;

                reset(&ramified_roots);
                mpz_set_ui(tmp_mpz, 0);
                append(&ramified_roots, tmp_mpz);

                mpz_set_ui(tmp_mpz, p);

                for (size_t j = 2 ; j < upto ; j++)
                {
                    dyn_array new;
                    init(&new);

                    for (size_t r = 0 ; r < ramified_roots.len ; r++)
                    {
                        mpz_mul_ui(tmp_mpz2, tmp_mpz, p);

                        evaluate_mod_p_mpz(f, tmp_mpz3, ramified_roots.start[r], tmp_mpz2);

                        if (!mpz_cmp_ui(tmp_mpz3, 0))
                        {
                            mpz_set(tmp_mpz2, ramified_roots.start[r]);
                            for (size_t k = 0 ; k < p ; k++)
                            {
                                append(&new, tmp_mpz2);
                                mpz_add(tmp_mpz2, tmp_mpz2, tmp_mpz);
                            }

                            mpf_set_z(tmp_mpf, tmp_mpz);
                            mpf_set_d(tmp_mpf2, 1.0);
                            mpf_div(tmp_mpf, tmp_mpf2, tmp_mpf);

                            Ep += mpf_get_d(tmp_mpf);
                            Fp += (2*j - 1)*mpf_get_d(tmp_mpf);
                        }
                    }

                    mpz_mul_ui(tmp_mpz, tmp_mpz, p);
                    copy_dyn(&ramified_roots, &new);

                    free_dyn_array(&new);
                }

                mpz_set_ui(tmp_mpz2, p);
                mpz_pow_ui(tmp_mpz2, tmp_mpz2, upto - 2);
                mpz_mul_ui(tmp_mpz2, tmp_mpz2, p - 1);

                mpf_set_z(tmp_mpf, tmp_mpz2);
                mpf_set_ui(tmp_mpf2, ramified_roots.len);
                mpf_div(tmp_mpf, tmp_mpf2, tmp_mpf);
                Ep += mpf_get_d(tmp_mpf);

                mpf_set_d(tmp_mpf, (double)upto*upto + 2*upto/(double)(p-1) + (p+1)/(double)((p-1)*(p-1)));
                mpf_set_z(tmp_mpf2, tmp_mpz);
                mpf_div(tmp_mpf, tmp_mpf, tmp_mpf2);
                mpf_mul_ui(tmp_mpf, tmp_mpf, ramified_roots.len);

                Fp += mpf_get_d(tmp_mpf);
            }

            free_polynomial(&frev);
            free_polynomial(&frev_derivative);
        }

        double tmpE = Ep*p/(p+1);
        double tmpF = Fp*p/(p+1);

        E += log_p * tmpE;
        F += (tmpF - tmpE * tmpE) * log_p * log_p;

        free(roots.start);
        free_dyn_array(&ramified_roots);
    }

    free_polynomial(&f_derivative);

    double k = 2*E - F/2;
    double lbd = F/2 - E;

    res[0] = E - baseline_term;
    res[1] = k;
    res[2] = lbd;
    res[3] = baseline_term;

    mpf_clears(tmp_mpf, tmp_mpf2, NULL);
    mpz_clears(tmp_mpz, tmp_mpz2, tmp_mpz3, NULL);
}

void compute_m_mu(
    mpz_t res,
    const mpz_t m0,
    const size_t l,
    const unsigned long d,
    const mpz_t components[l][d],
    const unsigned long * restrict vec
)
{
    mpz_set(res, m0);

    for (size_t i = 0 ; i < l ; i++)
    {
        mpz_add(res, res, components[i][vec[i]]);
    }
}

void compute_e(
    const mpz_t m0,
    const unsigned long nb_roots,
    const unsigned long d,
    const mpz_t roots_used[nb_roots][d],
    const mpz_t prod,
    const mpz_t a_d,
    const mpz_t n,
    mpz_t e[nb_roots][d]
)
{
    mpz_t tmp_m;
    mpz_init(tmp_m);

    unsigned long * vec = calloc(nb_roots, sizeof(unsigned long));

    polynomial_mpz tmp_poly;
    init_poly(&tmp_poly);

    compute_m_mu(tmp_m, m0, nb_roots, d, roots_used, vec);
    
    build_poly_coeffs(&tmp_poly, m0, tmp_m, a_d, n, d);

    mpz_t base;
    mpz_init_set(base, tmp_poly.coeffs[1]);
    mpz_mod(base, base, prod);

    for (size_t i = 0 ; i < nb_roots ; i++)
    {
        for (size_t j = 0 ; j < d ; j++)
        {
            if (!i)
            {
                for (size_t k = 0 ; k < nb_roots ; k++) vec[k] = 0;
                vec[0] = j;

                compute_m_mu(tmp_m, m0, nb_roots, d, roots_used, vec);
                build_poly_coeffs(&tmp_poly, m0, tmp_m, a_d, n, d);

                mpz_mod(e[i][j], tmp_poly.coeffs[1], prod);
            }
            else if (j)
            {
                for (size_t k = 0 ; k < nb_roots ; k++) vec[k] = 0;
                vec[i] = j;

                compute_m_mu(tmp_m, m0, nb_roots, d, roots_used, vec);
                build_poly_coeffs(&tmp_poly, m0, tmp_m, a_d, n, d);

                mpz_sub(e[i][j], tmp_poly.coeffs[1], base);
                mpz_mod(e[i][j], e[i][j], prod);
            }
            else
            {
                mpz_set_ui(e[i][0], 0);
            }
        }
    }


    mpz_clears(tmp_m, base, NULL);
    free_polynomial(&tmp_poly);

    free(vec);
}

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
)
{
    mpz_t tmp_mpz, prod2;
    mpz_inits(tmp_mpz, prod2, NULL);

    mpz_mul(prod2, prod, prod);

    mpf_t tmp_mpf, tmp_mpf2, tmp_mpf3;
    mpf_inits(tmp_mpf, tmp_mpf2, tmp_mpf3, NULL);

    mpz_pow_ui(tmp_mpz, m0, d);
    mpz_mul(tmp_mpz, tmp_mpz, a_d);
    mpz_sub(tmp_mpz, n, tmp_mpz);
    mpf_set_z(tmp_mpf, tmp_mpz);

    mpz_pow_ui(tmp_mpz, m0, d - 1);
    mpz_mul(tmp_mpz, tmp_mpz, prod);
    mpz_mul(tmp_mpz, tmp_mpz, prod);
    mpf_set_z(tmp_mpf2, tmp_mpz);

    mpf_div(f0, tmp_mpf, tmp_mpf2);

    for (size_t i = 0 ; i < nb_roots ; i++)
    {
        for (size_t j = 0 ; j < d ; j++)
        {
            mpz_mul_ui(tmp_mpz, a_d, d);
            mpz_mul(tmp_mpz, tmp_mpz, roots_used[i][j]);
            mpf_set_z(tmp_mpf, tmp_mpz);

            mpf_set_z(tmp_mpf2, prod2);

            mpf_div(tmp_mpf, tmp_mpf, tmp_mpf2);


            mpf_set_z(tmp_mpf2, e[i][j]);
            mpf_set_z(tmp_mpf3, prod);
            mpf_div(tmp_mpf2, tmp_mpf2, tmp_mpf3);

            mpf_add(tmp_mpf, tmp_mpf, tmp_mpf2);
            mpf_neg(f[i][j], tmp_mpf);
        }
    }

    mpz_clears(tmp_mpz, prod2, NULL);

    mpf_clears(tmp_mpf, tmp_mpf2, tmp_mpf3, NULL);
}

void create_first_array(
    const unsigned long nb_roots,
    const unsigned long d,
    const mpf_t f0,
    const mpf_t f[nb_roots][d],
    const unsigned long nb_rows,
    const unsigned long nb_cols,
    mpf_t * restrict array1_u_values,
    unsigned long array1_indices[nb_rows][nb_cols],
)
{
    unsigned long * vec = calloc(nb_cols, sizeof(unsigned long));

    size_t i = 0;

    mpf_t U, tmp_mpf;
    mpf_inits(U, tmp_mpf, NULL);

    while (vec[nb_cols  - 1] < d)
    {
        mpf_set(U, f0);

        for (size_t j = 0 ; j < nb_cols ; j++)
        {
            mpf_add(U, U, f[j][vec[j]]);
        }
        mpf_floor(tmp_mpf, U);
        mpf_sub(U, U, tmp_mpf);

        if (!i || mpf_cmp(U, array1_u_values[i]) > 0)
        {
            mpf_set(array1_u_values[i], U);
            for (size_t j = 0 ; j < nb_cols; j++) array1_indices[i][j] = vec[j];
        }
        else
        {
            unsigned long tmp_a = 0;
            unsigned long tmp_b = i - 1;
            unsigned long tmp_c = (tmp_a + tmp_b)>>1;

            while (tmp_a <= tmp_b)
            {
                if (mpf_cmp(U, array1_u_values[tmp_c]) < 0) tmp_b = tmp_c - 1;
                else tmp_a = tmp_c + 1;
                tmp_c = (tmp_a + tmp_b)>>1;
            }
            for (size_t j = i-1 ; j > tmp_a ; j++)
            {
                mpf_set(array1_u_values[j], array1_u_values[j-1]);
                for (size_t k = 0 ; k < nb_cols ; k++)
                {
                    array1_indices[j][k] = array1_indices[j-1][k];
                }
            }
            mpf_set(array1_u_values[tmp_a], U);
            for (size_t j = 0 ; j < nb_cols ; j++) array1_indices[tmp_a][j] = vec[j];
        }

        vec[0]++;
        for (size_t j = 0 ; j < nb_cols - 1 ; j++)
        {
            if (vec[j] == d)
            {
                vec[j] = 0;
                vec[j+1]++;
            }
            else break;
        }

        i++;
    }

    mpf_clears(U, tmp_mpf, NULL);

    free(vec);
}

void create_second_array(
    const unsigned long nb_roots,
    const unsigned long vec_len,
    const unsigned long d,
    const unsigned long nb_rows,
    const unsigned long nb_cols,
    const mpf_t f[nb_roots][d],
    mpf_t * restrict array2_u_values,
    unsigned long array2_indices[nb_rows][nb_cols]
)
{
    unsigned long * vec = calloc(nb_cols, sizeof(unsigned long));

    size_t i = 0;

    mpf_t U, tmp_mpf;
    mpf_inits(U, tmp_mpf, NULL);

    while (vec[nb_cols  - 1] < d)
    {
        mpf_set_ui(U, 0);

        for (size_t j = 0 ; j < nb_cols ; j++)
        {
            mpf_add(U, U, f[j + vec_len][vec[j]]);
        }
        mpf_neg(U, U);

        mpf_floor(tmp_mpf, U);
        mpf_sub(U, U, tmp_mpf);

        if (!i || mpf_cmp(U, array2_u_values[i]) > 0)
        {
            mpf_set(array2_u_values[i], U);
            for (size_t j = 0 ; j < nb_cols; j++) array2_indices[i][j] = vec[j];
        }
        else
        {
            unsigned long tmp_a = 0;
            unsigned long tmp_b = i - 1;
            unsigned long tmp_c = (tmp_a + tmp_b)>>1;

            while (tmp_a <= tmp_b)
            {
                if (mpf_cmp(U, array2_u_values[tmp_c]) < 0) tmp_b = tmp_c - 1;
                else tmp_a = tmp_c + 1;
                tmp_c = (tmp_a + tmp_b)>>1;
            }
            for (size_t j = i-1 ; j > tmp_a ; j++)
            {
                mpf_set(array2_u_values[j], array2_u_values[j-1]);
                for (size_t k = 0 ; k < nb_cols ; k++)
                {
                    array2_indices[j][k] = array2_indices[j-1][k];
                }
            }
            mpf_set(array2_u_values[tmp_a], U);
            for (size_t j = 0 ; j < nb_cols ; j++) array2_indices[tmp_a][j] = vec[j];
        }

        vec[0]++;
        for (size_t j = 0 ; j < nb_cols - 1 ; j++)
        {
            if (vec[j] == d)
            {
                vec[j] = 0;
                vec[j+1]++;
            }
            else break;
        }

        i++;
    }

    mpf_clears(U, tmp_mpf, NULL);

    free(vec);
}