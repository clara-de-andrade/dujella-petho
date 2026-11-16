#ifndef POLY_RING_H
#define POLY_RING_H

/*****************************************************************************
 * poly_ring.h implements methods for the poly_t type representing the usual *
 * operations on a ring of polynomials: negation, addition and subtraction,  *
 * scalar multiplication, polynomial multiplication (Karatsuba's algorithm), *
 * exponentiation by nonnegative integers, "shifting" exponents (equivalent  *
 * to multiplying a polynomial P(x) by some power x**p of x by a nonnegative *
 * integer p), evaluation at numerical values, polynomial composition, and   *
 * symbolic differentiation.                                                 *
 *****************************************************************************/

#include "poly_type.h"

#include <gmp.h>

/* ************************************************************************* */

// Sets `rop` to the negation of `op`.
void poly_neg(poly_t rop, poly_t op);
// Sets `rop` to the sum of `op1` and `op2`.
void poly_add(poly_t rop, poly_t op1, poly_t op2);
// Sets `rop` to the difference of `op1` and `op2`.
void poly_sub(poly_t rop, poly_t op1, poly_t op2);

/* ************************************************************************* */

// Sets `rop` to the product of `op1` and `op2`. Computed for small degrees
// of `op1` and `op2` by the Cauchy product formula for formal power series,
// and by Karatsuba's algorithm for larger degrees.
void poly_mul(poly_t rop, poly_t op1, poly_t op2);

// Sets `rop` to the scalar product of `op1` by `op2`.
void poly_mul_fr(poly_t rop, poly_t op1, mpfr_t op2);
void poly_mul_z(poly_t rop, poly_t op1, mpz_t op2);
void poly_mul_ui(poly_t rop, poly_t op1, unsigned long op2);
void poly_mul_si(poly_t rop, poly_t op1, long op2);

/* ************************************************************************* */

// Equivalent to setting `rop(x)` to the product of `op1(x)` by `x ** op2`.
void poly_shift(poly_t rop, poly_t op1, size_t op2);

// Sets `rop` to the power of `op1` by a non-negative integer `op2`.
void poly_pow_z(poly_t rop, poly_t op1, mpz_t op2);
// Sets `rop` to the power of `op1` by a non-negative integer `op2`.
// Faster than `poly_pow_z` for `unsigned long`s.
void poly_pow_ui(poly_t rop, poly_t op1, unsigned long op2);

/* ************************************************************************* */

// Sets `rop` to the composition of `op1` and `op2`.
void poly_eval(poly_t rop, poly_t op1, poly_t op2);

// Evaluates `op1(x)` at `x = op2` using Horner' method and sets to `rop.
void poly_eval_fr(mpfr_t rop, poly_t op1, mpfr_t op2);
void poly_eval_z(mpfr_t rop, poly_t op1, mpz_t op2);
void poly_eval_ui(mpfr_t rop, poly_t op1, unsigned long op2);
void poly_eval_si(mpfr_t rop, poly_t op1, long op2);

/* ************************************************************************* */

// Sets `rop` to derivative of `op`, computed by symbolic differentiation.
void poly_d(poly_t rop, poly_t op);

/* ************************************************************************* */

#endif