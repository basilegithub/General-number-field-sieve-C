#include <stdlib.h>

#include "dynamic_arrays.h"
#include "algebraic_base.h"

void algebraic_base_init(algebraic_base *b)
{
    b->start = NULL;
    b->end = NULL;
}

void algebraic_base_clear(algebraic_base *b)
{
    algebraic_base_prime *p = b->start;
    while (p != NULL) {
        algebraic_base_prime *next = p->next;

        free_dyn_array(&p->roots);

        free(p);

        p = next;
    }

    /* Reset the base */
    b->start = NULL;
    b->end = NULL;
}