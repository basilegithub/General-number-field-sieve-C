#ifndef SINGLE_LINKED_LIST_H
#define SINGLE_LINKED_LIST_H

typedef struct list_node
{
    unsigned long p;
    unsigned long r;
    struct list_node *next;
} list_node; // Dynamic algebraic base prime

typedef struct
{
    list_node *start;
    list_node *end;
} single_linked_list; // Dynamic algebraic base

// Functions declaration

void list_init(single_linked_list * restrict list);

void list_append(single_linked_list * restrict list, const unsigned long p, const unsigned long r);

void list_clear(single_linked_list * restrict list);

#endif // SINGLE_LINKED_LIST_H