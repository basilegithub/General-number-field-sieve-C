#include <gmp.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "algebraic_base.h"
#include "quadratic_characters.h"
#include "NFS_relations.h"

void mono_cpu_sieve(
    nfs_relations *relations,
    polynomial_mpz f_x,
    polynomial_mpz g_x,
    dyn_array_classic rat_base,
    algebraic_base alg_base,
    size_t nb_Algebraic_pairs,
    quadratic_character_base quad_base,
    size_t nb_Quadratic_characters,
    mpz_t leading_coeff,
    mpz_t prod_primes,
    mpz_t m0,
    mpz_t m1,
    mpz_t sieve_len,
    mpz_t const1,
    mpz_t const2,
    unsigned long *divide_leading,
    mpz_t *pow_div,
    size_t len_divide_leading,
    dyn_array_classic logs,
    FILE *logfile
)
{
    size_t required_relations = 3 + rat_base.len + nb_Algebraic_pairs + nb_Quadratic_characters;

    unsigned long offset = mpz_sizeinbase(const2, 2);

    log_msg(logfile, "Sieving started...");
    log_msg(logfile, "Need to collect at least %zu relations.", required_relations);

    unsigned long b = 1;

    mpz_t gcd_b_cd;
    mpz_init(gcd_b_cd);

    while (relations->len < required_relations)
    {
        // Update contribution of gcd of b and c_d

        mpz_gcd_ui(gcd_b_cd, leading_coeff, b);
        unsigned long new_offset = offset + mpz_sizeinbase(gcd_b_cd, 2);

        // Sieve

        // Verify smooth candidates

        // Append true smooths to the collected relations

        // Handle partial relations

        // Increment b

        b++;
    }
}