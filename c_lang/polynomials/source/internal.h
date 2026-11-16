#ifndef POLY_INTERNAL_H
#define POLY_INTERNAL_H

#include "polynomials.h"

#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>

/* ************************************************************************* */

// Largest degree until multiplication is computed by Karatsuba's algorithm.
#define POLY_MUL_CUTOFF 16

/* ************************************************************************* */

// Evaluates to maximum of two numbers, a and b.
#define max(a, b) ((a) >= (b) ? (a) : (b))
// Evaluates to minimum of two numbers, a and b.
#define min(a, b) ((a) <= (b) ? (a) : (b))

/* ************************************************************************* */

// Initialises guard for poly_t aliasing (when an operand is passed as `rop`).
#define POLY_GUARD_INIT(ptr, scratch, rop, cond) \
    int _aliased_##scratch = (cond);             \
    poly_t scratch;                              \
    struct _poly_t* ptr = (rop);                 \
    if (_aliased_##scratch)                      \
    {                                            \
        poly_init(scratch);                      \
        ptr = scratch;                           \
    }

// Clears guard for poly_t aliasing (when an operand is passed as `rop`).
#define POLY_GUARD_CLEAR(rop, scratch)           \
    if (_aliased_##scratch)                      \
    {                                            \
        poly_set((rop), scratch);                \
        poly_clear(scratch);                     \
    }

/* ************************************************************************* */

// Initialises guard for mpfr_t aliasing
#define MPFR_GUARD_INIT(ptr, scratch, rop, cond)    \
    int _aliased_##scratch = (cond);                \
    mpfr_t scratch;                                 \
    mpfr_ptr ptr = (rop);                           \
    if (_aliased_##scratch)                         \
    {                                               \
        mpfr_init2(scratch, mpfr_get_prec(rop));    \
        mpfr_set(scratch, (rop), MPFR_RNDN);        \
        ptr = scratch;                              \
    }

// Clears guard for mpfr_t aliasing
#define MPFR_GUARD_CLEAR(scratch)   \
    if (_aliased_##scratch)         \
        mpfr_clear(scratch);

/* ************************************************************************* */

// Checks whether `ptr` points within array `base` of given length `len`.
#define ptr_range(ptr, base, len)                           \
    ((base) != NULL &&                                      \
     (uintptr_t) (ptr) >= (uintptr_t) (base) &&             \
     (uintptr_t) (ptr) <  (uintptr_t) ((base) + (len)))

// Initialises guard for poly_t scalar aliasing
#define POLY_SCALAR_GUARD_INIT(ptr, scratch, rop, op) \
    MPFR_GUARD_INIT(ptr, scratch, op, ptr_range(op, rop->coeffs, rop->alloc))

// Clears guard for poly_t scalar aliasing
#define POLY_SCALAR_GUARD_CLEAR(scratch) \
    MPFR_GUARD_CLEAR(scratch)

/* ************************************************************************* */

// Reallocates `rop->coeffs` to either twice as many or `alloc` elements.
void poly_realloc(poly_t rop, size_t alloc);

// Sets `rop->deg` to `deg`. If needed, reallocates or truncates `rop->coeff`.
void poly_set_deg(poly_t rop, size_t deg);

// Sets `rop->deg` to largest `index` such that `rop->coeff[index] != 0`.
void poly_normalise(poly_t rop);

/* ************************************************************************* */

#endif
