#ifndef DUJELLA_PETHO_H
#define DUJELLA_PETHO_H

/*****************************************************************************
 * dujella_petho.h contains a C type, dp_t, encoding a structure with the    *
 * parameters that appear in the Dujella-Petho Theorem, whose statement is:  *
 *                                                                           *
 * Let M >= be an integer and gamma an irrational number; let p/q be any     *
 * convergent of the continued fraction of gamma such that q > 6*M; and      *
 * finally, let mu be an irrational number and put                           *
 *                                                                           *
 *   eps = ||mu * q|| - M * ||gamma * q||,                                   *
 *                                                                           *
 * where ||-|| denotes the distance to the nearest integer. If eps > 0 then  *
 * there does not exist any integer solutions to the inequality              *
 *                                                                           *
 *   0 < |m * gamma - k + mu| < A * B**-n                                    *
 *                                                                           *
 * with m <= M and n >= N, where N = log(A * q / eps) / log(B).              *
 *                                                                           *
 * dujella_petho.h also implements methods to initialise and clear each of   *
 * these parameters, and a function that, given gamma, mu, A, B and M,       *
 * computes the least q such that q > 6*M and eps > 0, as well as eps and    *
 * the upper bound N.                                                        *
 *****************************************************************************/

 /* TODO:
  * [ ] Document dujella_petho.c.
  * [ ] Write computation of q and eps separately.
  */

#include <gmp.h>
#include <mpfr.h>

/* ************************************************************************* */

// C struct of Dujella-Petho parameters
struct _dp_t {
    mpfr_prec_t prec;
    
    mpfr_t gamma;
    mpfr_t mu;
    mpfr_t A;
    mpfr_t B;
    mpfr_t eps;

    mpz_t M;
    mpz_t q;
    mpz_t N;
};

// C type of Dujella-Petho parameters
typedef struct _dp_t dp_t[1];

/* ************************************************************************* */

// Initializes parameters of `rop` with precision `prec`.
void dp_init(dp_t rop, mpfr_prec_t prec);

// Clears parameters of `rop`.
void dp_clear(dp_t rop);

/* ************************************************************************* */

// Computes dp->q, dp->eps and dp->N from remaining parameters.
void dujella_petho(dp_t dp);

/* ************************************************************************* */

#endif