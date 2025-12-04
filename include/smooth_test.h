#ifndef SMOOTH_TEST_H
#define SMOOTH_TEST_H

#include <gmp.h>

#include "dynamic_arrays.h"
#include "single_linked_list.h"
#include "polynomial_structures.h"
#include "NFS_relations.h"

void pollard_rho(const mpz_t m, mpz_t p1, mpz_t p2, gmp_randstate_t state);
void naive_smooth(nfs_relations *smooth_candidates, dyn_array_classic primes, mpz_t limit, mpz_t limit_2, gmp_randstate_t state);


#endif // SMOOTH_TEST_H