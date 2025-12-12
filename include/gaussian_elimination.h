#ifndef GAUSSIAN_ELIMINATION_H
#define GAUSSIAN_ELIMINATION_H

#include <stdbool.h>

void gaussian_elimination(mpz_t * restrict dense_matrix, mpz_t * restrict res, const unsigned long relations_len, const unsigned long base_size);
bool row_is_zero(const mpz_t * restrict dense_matrix, const size_t row_index);

#endif // GAUSSIAN_ELIMINATION_H