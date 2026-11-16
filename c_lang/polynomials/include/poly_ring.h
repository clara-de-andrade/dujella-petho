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

// Tests if `op` is zero, where `op` is a scalar.
#define mpfr_zero(op) (mpfr_zero_p(op))
#define mpz_zero(op)  (mpz_sgn(op) == 0)
#define ui_zero(op)   (op == 0)
#define si_zero(op)   (op == 0)

// Tests if `op` is zero, where `op` is a polynomial.
#define poly_zero(op) ((op)->deg == 0 && mpfr_zero((op)->coeffs[0]))


// Tests if `op` is one, where `op` is a scalar.
#define mpfr_unit(op) (!mpfr_nan_p(op) && mpfr_cmp_ui(op, 1) == 0)
#define mpz_unit(op)  (mpz_cmp_ui(op, 1) == 0)
#define ui_unit(op)   (op == 1)
#define si_unit(op)   (op == 1)

// Tests if `op` is one, where `op` is a polynomial.
#define poly_unit(op) ((op)->deg == 0 && mpfr_unit((op)->coeffs[0]))

/* ************************************************************************* */

// Sets `rop` to the negation of `op`.
void poly_neg(poly_t rop, poly_t op);
// Sets `rop` to the sum of `op1` and `op2`.
void poly_add(poly_t rop, poly_t op1, poly_t op2);
// Sets `rop` to the difference of `op1` and `op2`.
void poly_sub(poly_t rop, poly_t op1, poly_t op2);

/* ************************************************************************* */

// Sets `rop` to the product of `op1` and `op2` modulo `x**d`.
void poly_mullow(poly_t rop, poly_t op1, poly_t op2, size_t d);

// Sets `rop` to the product of `op1` and `op2`.
#define poly_mul(rop, op1, op2) \
    poly_mullow(rop, op1, op2, op1->deg + op2->deg + 1)

// Sets `rop(x)` to the product of `op1(x)` by `x ** op2`.
void poly_mul_exp(poly_t rop, poly_t op1, size_t op2);

// Sets `rop` to the scalar product of `op1` by `op2`.
void poly_mul_fr(poly_t rop, poly_t op1, mpfr_t op2);
void poly_mul_z(poly_t rop, poly_t op1, mpz_t op2);
void poly_mul_ui(poly_t rop, poly_t op1, unsigned long op2);
void poly_mul_si(poly_t rop, poly_t op1, long op2);

/* ************************************************************************* */

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