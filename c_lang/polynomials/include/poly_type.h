#ifndef POLY_TYPE_H
#define POLY_TYPE_H

/*****************************************************************************
 * poly_type.h contains a C Type named poly_t, implementing polynomials with *
 * MPFR floating numbers as coefficients as well as methods for operating on *
 * these objects as a class; namely, methods for initialising, dedclaring,   *
 * assigning and clearing polynomials in memory, setter and getter methods   *
 * for public attributes (coefficients and degree), and printing.            *
 *****************************************************************************/

#include <stddef.h>
#include <mpfr.h>

/*****************************************************************************
 * Polynomials are represented by "sequences" (i.e. arrays) of coefficients, *
 * which themselves are MPFR floating point numbers. This array is allocated *
 * when the polynomial is initialised, and the coefficients are set to zero  *
 * by default. More concretely, a poly_t is a pointer to a struct _poly_t,   *
 * which contains an allocated array of coefficients, as well as the length  *
 * of this array and the degree of the polynomial represented by the poly_t. *
 *                                                                           *
 * Note that we define poly_t as a one-element array of struct _poly_t's, as *
 * in this case it doesn't need allocation and can be handled by the stack.  *
 * However, we effectively treat it as a pointer to a single struct _poly_t. *
 *****************************************************************************/

// C struct of polynomials
struct _poly_t {
    mpfr_t* coeffs;     // Array of coefficients of polynomial
    size_t alloc;       // Number of allocated coefficients
    size_t deg;         // Degree of polynomial
};

// C type of polynomials
typedef struct _poly_t poly_t[1];

/* ************************************************************************* */

// Initializes polynomial `rop` to zero polynomial.
void poly_init(poly_t rop);
// Initializes each argument with `poly_init`.
void poly_inits(poly_t rop, ...);

// Clears polynomial `rop`.
void poly_clear(poly_t rop);
// Clears each argument with `poly_clear`.
void poly_clears(poly_t rop, ...);

/* ************************************************************************* */

// Declares polynomial given array of coefficients and degree.
void poly_new_fr(poly_t rop, mpfr_t* coeffs, size_t deg);
void poly_new_z(poly_t rop, mpz_t* coeffs, size_t deg);
void poly_new_si(poly_t rop, const long* coeffs, size_t deg);
void poly_new_ui(poly_t rop, const unsigned long* coeffs, size_t deg);

// Initialises and declares polynomial given array of coefficients and degree.
#define poly_init_new_fr(rop, coeffs, deg)  \
    do {                                    \
        poly_init(rop);                     \
        poly_new_fr(rop, coeffs, deg);      \
    } while (0)
#define poly_init_new_z(rop, coeffs, deg)  \
    do {                                    \
        poly_init(rop);                     \
        poly_new_z(rop, coeffs, deg);      \
    } while (0)
#define poly_init_new_si(rop, coeffs, deg)  \
    do {                                    \
        poly_init(rop);                     \
        poly_new_si(rop, coeffs, deg);      \
    } while (0)
#define poly_init_new_ui(rop, coeffs, deg)  \
    do {                                    \
        poly_init(rop);                     \
        poly_new_ui(rop, coeffs, deg);      \
    } while (0)

/* ************************************************************************* */

// Sets polynomial `rop` to a copy of polynomial `op`.
void poly_set(poly_t rop, poly_t op);

// Initialises `rop` and sets to a copy of `op`.
#define poly_init_set(rop, op)      \
    do {                            \
        poly_init(rop);             \
        poly_set(rop, op);          \
    } while (0)


// Sets `rop` to constant polynomial `rop(x) = op`, where `op` is a scalar.
void poly_set_fr(poly_t rop, mpfr_t op);
void poly_set_z(poly_t rop, mpz_t op);
void poly_set_ui(poly_t rop, unsigned long op);
void poly_set_si(poly_t rop, long op);

// Initialises `rop` and sets to a constant polynomial `rop(x) = op`.
#define poly_init_set_fr(rop, op)   \
    do {                            \
        poly_init(rop);             \
        poly_set_fr(rop, op);       \
    } while (0)
#define poly_init_set_z(rop, op)    \
    do {                            \
        poly_init(rop);             \
        poly_set_z(rop, op);        \
    } while (0)
#define poly_init_set_ui(rop, op)   \
    do {                            \
        poly_init(rop);             \
        poly_set_ui(rop, op);       \
    } while (0)
#define poly_init_set_si(rop, op)   \
    do {                            \
        poly_init(rop);             \
        poly_set_si(rop, op);       \
    } while (0)

/* ************************************************************************* */

// Returns `rop->deg`.
size_t poly_get_deg(poly_t rop);
// Sets `rop` to `op->coeffs[index]`, setting to `0` if `index > op->deg`.
void poly_get_coeff(mpfr_t rop, poly_t op, size_t index);

// Sets `rop->coeffs[index]` to `op`, reallocating `rop` if `index > op->deg`.
void poly_set_coeff_fr(poly_t rop, mpfr_t op, size_t index);
void poly_set_coeff_z(poly_t rop, mpz_t op, size_t index);
void poly_set_coeff_ui(poly_t rop, unsigned long op, size_t index);
void poly_set_coeff_si(poly_t rop, long op, size_t index);

/* ************************************************************************* */

// Returns string representation of polynomial `op`.
char* poly_repr(poly_t op);

/* ************************************************************************* */

#endif