#ifndef INIT_FUNCTIONS_H
#define INIT_FUNCTIONS_H

#include <gmp.h>

unsigned int compute_degree(mpz_t n, mpf_t ln2, mpf_t e);
void initialize_params(gmp_randstate_t state, mpf_t ln10, mpf_t ln2, mpf_t e);
void compute_smooth_bound(mpz_t n, mpz_t smooth_bound, mpf_t ln2, mpf_t e);

#endif // INIT_FUNCTIONS_H