#ifndef __POLYNOM_H__
#define __POLYNOM_H__

#include <stdint.h>
#include <mpfr.h>
#include <gmp.h>

/* ************************************************************************* */

// C Type of polynomials with integer coefficients
typedef struct _poly_t poly_t;

// Creates a new polynomial coeffs[0] + ... + coeffs[deg] * x**(deg)
poly_t* poly_new(size_t deg, mpz_t* coeffs);
// Frees self
void poly_free(poly_t* self);

/* ************************************************************************* */

// Getter for deg
size_t poly_get_deg(poly_t* self);
// Getter for coeffs
mpz_t* poly_get_coeffs(poly_t* self);
// Getter for coeffs[key]
void poly_get_coeff(poly_t* self, size_t key, mpz_t op);

// Setter for deg
void poly_set_deg(poly_t* self, size_t value);
// Setter for coeffs
void poly_set_coeffs(poly_t* self, mpz_t* value);
// Setter for coeffs[key]
void poly_set_coeff(poly_t* self, size_t key, mpz_t value);

/* ************************************************************************* */

// String representation of self
char* poly_str(poly_t* self);
// Number of coefficients in self
size_t poly_len(poly_t* self);

// Evaluate self(x) at x = other
int64_t poly_eval_si(poly_t* self, int64_t other);
uint64_t poly_eval_ui(poly_t* self, uint64_t other);
void poly_eval_z(poly_t* self, mpz_t other, mpz_t op);
void poly_eval_fr(poly_t* self, mpfr_t other, mpfr_t op);
poly_t* polt_eval(poly_t* self, poly_t* other);

/* ************************************************************************* */

// Computes self(x) + other(x) by pointwise addition
poly_t* poly_add(poly_t* self, poly_t* other);
// Computes self(x) - other(x) by pointwise addition
poly_t* poly_sub(poly_t* self, poly_t* other);

// Computes other * self(x) by scalar multiplication
poly_t* poly_mul_si(poly_t* self, int64_t other);
poly_t* poly_mul_ui(poly_t* self, uint64_t other);
poly_t* poly_mul_z(poly_t* self, mpz other);
poly_t* poly_mul_fr(poly_t* self, mpfr_t other);
// Computes self(x) * other(x) using Karatsuba's algorithm
poly_t* poly_mul(poly_t* self, poly_t* other);

// Computes self(x) ** exp using exponentiation by squaring
poly_t* poly_pow_ui(poly_t* self, uint64_t exp);
poly_t* poly_pow_z(poly_t* self, mpz_t exp);

// Computes -self(x) = (-1) * self(x)
poly_t* poly_neg(poly_t* self);

/* ************************************************************************* */

// Differentiates of polynomial P(t) = a_0 + ... + a_n * t**n, where n >= 0.
poly_t* poly_diff(poly_t* P);

/* ************************************************************************* */

#endif
