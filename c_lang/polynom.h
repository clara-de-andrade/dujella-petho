#ifndef __POLYNOM_H__
#define __POLYNOM_H__

#include <stdint.h>
#include <mpfr.h>
#include <gmp.h>

/* ************************************************************************* */

// Creates polynomial P(x) = coeffs[0] + ... + coeffs[deg - 1] * x**(deg - 1).
poly_t* poly_new(size_t deg, mpz_t* coeffs);

// Frees poly_t pointer allocated during poly_init.
void poly_free(poly_t* poly);

/* ************************************************************************* */

size_t poly_get_deg(poly_t* poly);
mpz_t* poly_get_coeffs(poly_t* poly);
mpz_t poly_get_coeff(poly_t* poly, size_t n);

void poly_set_deg(poly_t* poly, size_t deg);
void poly_set_coeffs(poly_t* poly, mpz_t* coeffs);
void poly_set_coeff(poly_t* poly, size_t n, mpz_t coeff);

/* ************************************************************************* */

// Computes the sum of two polynomials A(t) = a_0 + ... + a_m * t**m and
// B = b_0 + ... + b_n * t**n, where m, n >= 0, by pointwise addition.
poly_t* poly_add(poly_t* A, poly_t* B);

// Computes the product of two polynomials A(t) = a_0 + ... + a_m * t**m and
// B = b_0 + ... + b_n * t**n, where m, n >= 0, using Karatsuba's algorithm.
poly_t* poly_mul(poly_t* A, poly_t* B);

/* ************************************************************************* */

// Differentiates of polynomial P(t) = a_0 + ... + a_n * t**n, where n >= 0.
poly_t* poly_diff(poly_t* P);

/* ************************************************************************* */

// Evaluates polynomial P(t) = a_0 + ... + a_n * t**n at real value t = x.
mpfr_t poly_eval(poly_t* P, mpfr_t x);

/* ************************************************************************* */

#endif
