#include <stdlib.h>

#include "dynamic_arrays.h"
#include "single_linked_list.h"

void list_init(single_linked_list *list)
{
    list->start = NULL;
    list->end = NULL;
}

void list_append(single_linked_list *list, unsigned long p, unsigned long r)
{
    list_node *next_node = malloc(sizeof(list_node));
    next_node->p = p;
    next_node->r = r;
    next_node->next = NULL;

    if (list->start == NULL)
    {
        list->start = next_node;
        list->end = next_node;
    }
    else
    {
        list->end->next = next_node;
        list->end = next_node;
    }
}

void list_clear(single_linked_list *list)
{
    list_node *p = list->start;
    while (p != NULL) {
        list_node *next = p->next;

        free(p);

        p = next;
    }

    list->start = NULL;
    list->end = NULL;
}