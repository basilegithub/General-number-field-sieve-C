#include <stdlib.h>

#include "dynamic_arrays.h"
#include "single_linked_list.h"
#include "NFS_relations.h"

void init_relations(nfs_relations *relations)
{
    relations->rels = calloc(1, sizeof(nfs_relation));
    relations->len = 0;
    relations->size = 1;
}

void init_new_relation(nfs_relations *relations)
{
    if (relations->len == relations->size)
    {
        relations->size <<= 1;
        relations->rels = realloc(relations->rels, relations->size * sizeof(nfs_relation));
    }

    poly_init(&relations->rels[relations->len].poly_g);
    poly_init(&relations->rels[relations->len].poly_f);
    mpz_init(relations->rels[relations->len].rational_norm);
    dyn_array_init(&relations->rels[relations->len].rational_large_primes);
    mpz_init(relations->rels[relations->len].algebraic_norm);
    list_init(relations->rels[relations->len].algebraic_large_primes);
    relations->rels[relations->len].nb_relations = 0;

    relations->len++;
}

void clear_relations(nfs_relations *relations)
{
    for (size_t i = 0 ; i < relations->len ; i++)
    {
        free_polynomial(&relations->rels[i].poly_g);

        free_polynomial(&relations->rels[i].poly_f);

        mpz_clear(relations->rels[i].rational_norm);

        free(relations->rels[i].rational_large_primes.start);
        relations->rels[i].rational_large_primes.start = NULL;
        relations->rels[i].rational_large_primes.len = 0;
        relations->rels[i].rational_large_primes.size = 0;

        mpz_clear(relations->rels[i].algebraic_norm);

        list_clear(&relations->rels[i].algebraic_large_primes);
    }
    
    free(relations->rels);
    relations->rels = NULL;

    relations->len = 0;
    relations->size = 0;
}