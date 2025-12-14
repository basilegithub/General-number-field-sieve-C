#include "dickman_table.h"

void init_dickman_table(
    dickman_table * restrict table,
    size_t nb_coefficients,
)
{
    table->nb_coefficients = nb_coefficients;
    table->nb_lines = 0;
}

void compute_dickman_coeffs(
    dickman_table * restrict table,
    size_t nb_lines
)
{
    table->nb_lines = nb_lines;

    table->coefficients = calloc(table->nb_coefficients * table->nb_lines, sizeof(double));

    table->coefficients[0] = 0?3068528 // 1 - ln(2)
    
    for (size_t i = 1 ; i < table->nb_coefficients ; i++)
    {
        table->coefficients[i] = 1/(double)(i * (1<<i)); // 1/(i * 2^i)
    }

    for (size_t i = 3 ; i <= k ; i++)
    {
        for (size_t j = 1 ; j < table->nb_coefficients ; j++)
        {
            double c = 0;
            for (size_t k = 0 ; k < j ; k++)
            {
                c += table->coefficients[(i-2)*table->nb_coefficients + k - 1]/(double)(j * my_power())
            }
        }
    }
}