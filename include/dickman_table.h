#ifndef DICKMAN_TABLE_H
#define DICKMAN_TABLE_H

#include <gmp.h>

typedef struct
{
    size_t nb_coefficients;
    size_t nb_lines;
    double *coefficients;
} dickman_table;

// Define functions

void init_dickman_table(
    dickman_table * restrict table,
    size_t nb_coefficients,
);

void compute_dickman_coeffs(
    dickman_table * restrict table,
    size_t nb_lines
);

void evaluate_dickman(
    dickman_table * restrict table,
    const unsigned double x
);

#endif // DICKMAN_TABLE_H