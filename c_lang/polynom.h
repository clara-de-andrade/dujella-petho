#ifndef __POLYNOM_H__
#define __POLYNOM_H__

#include <stdint.h>
#include <mpfr.h>
#include <gmp.h>

/* ************************************************************************* */

// C Type of polynomials with integer coefficients
typedef struct {
    size_t deg;
    mpz_t* coeffs;
} poly_t;

/* ************************************************************************* */

// Creates polynomial P(x) = coeffs[0] + ... + coeffs[deg - 1] * x**(deg - 1).
poly_t* poly_init(size_t deg, mpz_t* coeffs)

// Frees poly_t pointer allocated during poly_init.
void poly_clear(poly_t* P);

// Computes the sum of two polynomials A(t) = a_0 + ... + a_m * t**m and
// B = b_0 + ... + b_n * t**n, where m, n >= 0, by pointwise addition.
poly_t* poly_add(poly_t* A, poly_t* B);

// Computes the product of two polynomials A(t) = a_0 + ... + a_m * t**m and
// B = b_0 + ... + b_n * t**n, where m, n >= 0, using Karatsuba's algorithm.
poly_t* poly_mul(poly_t* A, poly_t* B);

// Differentiates of polynomial P(t) = a_0 + ... + a_n * t**n, where n >= 0.
poly_t* poly_diff(poly_t* P);

// Evaluates polynomial P(t) = a_0 + ... + a_n * t**n at real value t = x.
void poly_eval(mpfr_t op, poly_t* P, mpfr_t x);

/* ************************************************************************* */

#endif
