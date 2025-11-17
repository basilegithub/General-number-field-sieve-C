#include "quadratic_characters.h"

void algebraic_base_init(quadratic_character_base *b)
{
    b->start = NULL;
}

void algebraic_base_clear(quadratic_character_base *b)
{
    quadratic_character *p = b->start;
    while (p != NULL) {
        quadratic_character *next = p->next;

        free(p);

        p = next;
    }

    b->start = NULL;
}