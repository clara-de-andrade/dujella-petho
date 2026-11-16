#ifndef POLY_ROOT_H
#define POLY_ROOT_H

/*****************************************************************************
 * poly_root.h contains methods for finding approximations to the real roots *
 * of a polynomial in various cases: when a close approximation is known...  *
 *****************************************************************************/

#include "poly_type.h"

/* ************************************************************************* */

// TODO: more methods

#define POLY_ROOT_CONVERGENCE    0
#define POLY_ROOT_STATIONARY    -1
#define POLY_ROOT_MAXITER       +1

// Sets `rop` to root of `poly` closest to `init`. Returns 0 on convergence,
// positive if exceeded maximum number of iterations, and negative if either
// failed to convergence or found a stationary point. 
void poly_root(mpfr_t rop, poly_t poly, mpfr_t init);

/* ************************************************************************* */

#endif