#include <stdlib.h>

#include "dickman_table.h"
#include "utils.h"

void init_dickman_table(
    dickman_table * restrict table,
    const size_t nb_coefficients
)
{
    table->nb_coefficients = nb_coefficients;
    table->nb_lines = 0;
}

void compute_dickman_coeffs(
    dickman_table * restrict table,
    const size_t nb_lines
)
{
    table->nb_lines = nb_lines;

    table->coefficients = calloc(table->nb_coefficients * table->nb_lines, sizeof(double));

    table->coefficients[0] = 0.3068528; // 1 - ln(2)
    
    for (size_t i = 1 ; i < table->nb_coefficients ; i++)
    {
        table->coefficients[i] = 1/(double)(i * (1<<i)); // 1/(i * 2^i)
    }

    for (size_t i = 3 ; i <= nb_lines ; i++)
    {
        for (size_t j = 1 ; j < table->nb_coefficients ; j++)
        {
            double c = 0;
            for (size_t k = 0 ; k < j ; k++)
            {
                c += table->coefficients[(i-3)*table->nb_coefficients + k]/(double)(j * my_power(i, j - k));
            }
            table->coefficients[(i-2)*table->nb_coefficients + j] = c;
        }

        double c = 0;

        for (size_t j = 1 ; j < table->nb_coefficients ; j++)
        {
            c += table->coefficients[(i-2)*table->nb_coefficients + j]/(double)(j+1);
        }

        table->coefficients[(i-2)*table->nb_coefficients] = c/(i-1);
    }
}

double evaluate_dickman(
    dickman_table * restrict table,
    const double x
)
{
    unsigned long integer_part = (unsigned long) x;
    double delta = integer_part + 1 - x;

    if (integer_part > table->nb_lines)
    {
        free(table->coefficients);
        compute_dickman_coeffs(table, integer_part + 1);
    }

    double tmp = 1.0;
    double res = 0.0;

    for (size_t i = 0 ; i < table->nb_coefficients ; i++)
    {
        res += table->coefficients[table->nb_coefficients * (integer_part - 1) + i] * tmp;
        tmp *= delta;
    }

    return res;
}