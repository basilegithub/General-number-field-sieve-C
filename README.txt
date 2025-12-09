Thanks for looking at this project.

This is a new project aiming at concluding my work on factoring algorithms.


##### Introduction ######

This project is a C implementation of the General number field sieve to factor integers.

This a first working version, which contains the following features:
- Basic polynomial selection (d-th root) ;
- Naive sieving ;
- Naive and Batch smoothness test ;
- No large primes allowed ;
- Block Lanczos used for the linear algebra step (Gaussian elimination and Wiedemann algorithm implemented, but not available to choose for now) ;
- Lifting method to extract the algebraic square root.

##### Sources #####

Overall algorithm:
- "Prime numbers, a computational perspective" by Richard Crandall and Carl Pomerance (really good): https://link.springer.com/book/10.1007/0-387-28979-8
- "The Development of the Number Field Sieve" by many (really really good): https://link.springer.com/book/10.1007/BFb0091534
- "A beginner's guide to the general number field sieve" by Michael Case: https://www.cs.umd.edu/~gasarch/TOPICS/factoring/NFSmadeeasy.pdf

Kleinjung polynomials search algorithm:
- "On polynomial selection for the number field sieve" by THORSTEN KLEINJUNG: https://www.ams.org/journals/mcom/2006-75-256/S0025-5718-06-01870-9/S0025-5718-06-01870-9.pdf

Polynomial metric Ep_score:
- "A new ranking function for polynomial selection in the number field sieve" by Nicolas David and Paul Zimmerman: https://inria.hal.science/hal-02151093v4/document

Double large prime
- "Factoring with two large primes" by Arjen K. Lenstra and Mark S. Manasse: https://scispace.com/pdf/factoring-with-two-large-primes-1lk9719aco.pdf

Batch smoothness test:
- "HOW TO FIND SMOOTH PARTS OF INTEGERS" by DANIEL J. BERNSTEIN: https://cr.yp.to/factorization/smoothparts-20040510.pdf

Gaussian elimination:
- Took the gaussian elimination code from https://github.com/skollmann/PyFactorise/blob/master/factorise.py#L39

Block Lanczos:
- "A Block Lanczos Algorithm for Finding Dependencies over GF(2)" by Peter L. Montgomery: https://scispace.com/pdf/a-block-lanczos-algorithm-for-finding-dependencies-over-gf-2-ezdu2qt0pp.pdf
- "A modified block Lanczos algorithm with fewer vectors" by Emmanuel Thomé (really good): https://eprint.iacr.org/2016/329.pdf

Wiedemann algorithm:
- "SOLVING HOMOGENEOUS LINEAR EQUATIONSOVER GF(2) VIA BLOCK WIEDEMANN ALGORITHM" by Don Coppersmith: https://www.ams.org/journals/mcom/1994-62-205/S0025-5718-1994-1192970-7/S0025-5718-1994-1192970-7.pdf

Square root algorithms
- "Computing a square root for the number field sieve" by Jean-Marc Couveignes: https://www.math.u-bordeaux.fr/~jcouveig/publi/Cou94-2.pdf

##### running the algortihm #####

You can run the GNFS.exe in the build folder, it will ask you for the number you want to factor.

##### config parameters #####

None for now.

##### General discussion #####

Here are the next steps :

- Implement large primes handling ;

- Implement Kleinjung polynomial search algorithm ;

- Implement parallel polynomial search ;

- Implement parallel sieving ;

- Implement Couveignes algorithm for square root extraction ;

- Setup the config file ;

- Optimize heavily the code.