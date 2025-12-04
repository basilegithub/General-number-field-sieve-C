#include <gmp.h>

#include "dynamic_arrays.h"
#include "polynomial_structures.h"
#include "NFS_relations.h"
#include "utils.h"
#include "square_root.h"

void extract_solution(
    nfs_relations *relations,
    bool *kernel_vector,
    dyn_array_classic *rational_primes,
    polynomial_mpz f_x,
    polynomial_mpz f_prime_sq,
    mpz_t n,
    mpz_t m0,
    mpz_t m1,
    mpz_t f_prime_eval,
    unsigned long inert_prime,
    unsigned long max_a_size
)
{
    mpz_t f_norm;
    mpz_init(f_norm);

    for (size_t i = 0 ; i <= f_x.degree ; i++)
    {
        mpz_addmul(f_norm, f_x.coeffs[i], f_x.coeffs[i]);
    }
    mpz_sqrt(f_norm, f_norm);

    mpz_t fd;
    mpz_init_set_ui(fd, f_x.degree);
    mpz_pow_ui(fd, fd, 3);
    mpz_sqrt(fd, fd);

    mpz_t x, rational_square_root;
    mpz_init_set(x, f_prime_eval);
    mpz_init(rational_square_root);

    extract_rational_square_root(rational_square_root, relations, kernel_vector, n, rational_primes);
    
    gmp_printf("x = %Zd\n", rational_square_root);
    
    mpz_clears(f_norm, fd, NULL);
}