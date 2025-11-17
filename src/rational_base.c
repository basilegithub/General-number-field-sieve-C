#include <stdlib.h>

#include "dynamic_arrays.h"
#include "rational_base.h"

void rational_base_init(rational_base *b)
{
    b->start = NULL;
    b->end = NULL;
}

void rational_base_clear(rational_base *b)
{
    rational_base_prime *p = b->start;
    while (p != NULL) {
        rational_base_prime *next = p->next;

        free_dyn_array(&p->roots);

        free(p);

        p = next;
    }

    /* Reset the base */
    b->start = NULL;
    b->end = NULL;
}