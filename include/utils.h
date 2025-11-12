#ifndef UTILS_H
#define UTILS_H

#include <gmp.h>

void compute_e(mpf_t e);
void initialize_params(gmp_randstate_t state, mpf_t ln10, mpf_t ln2, mpf_t e);
void natural_log(mpf_t res, mpf_t x, mpf_t ln2, mpf_t e);

#endif // UTILS_H