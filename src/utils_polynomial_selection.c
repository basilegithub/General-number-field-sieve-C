#include <gmp.h>
#include <math.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "polynomial_functions.h"
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