#ifndef NFS_RELATIONS_H
#define NFS_RELATIONS_H

#include <gmp.h>

#include "dynamic_arrays.h"
#include "single_linked_list.h"
#include "polynomial_structures.h"

typedef struct
{
    polynomial_mpz poly_g;
    polynomial_mpz poly_f;
    mpz_t rational_norm;
    dyn_array_classic rational_large_primes;
    mpz_t algebraic_norm;
    single_linked_list algebraic_large_primes;
    unsigned long nb_relations; // For relations obtained from partials
    bool *divide_leading;

} nfs_relation; // Structure that holds all the required information of a nfs relation

typedef struct 
{
    nfs_relation *rels;
    size_t len;
    size_t size;
} nfs_relations;

// Functions declaration

void init_relations(nfs_relations * restrict relations);
void init_new_relation(nfs_relations * restrict relations, const size_t len_divide_leading);
void clear_relations(nfs_relations * restrict relations);


#endif // NFS_RELATIONS_H