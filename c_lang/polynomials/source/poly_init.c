#include "internal.h"

#include <stdarg.h>

/* ************************************************************************* */

// Initializes polynomial `rop` to zero polynomial, by allocating `rop->coeffs`
// and setting `rop->alloc` and `rop->deg` accordingly.
void poly_init(poly_t rop)
{
    // Allocate a one-element array of `mpfr_t`s
    rop->coeffs = malloc(sizeof(mpfr_t));

    // MEMORY ERROR
    if (rop->coeffs == (mpfr_t*) NULL)
        exit(EXIT_FAILURE);

    // Initialises constant term to zero by default
    mpfr_init_set_ui(rop->coeffs[0], 0, MPFR_RNDN);

    rop->alloc = 1;
    rop->deg = 0;
}

// Clears polynomial `rop` by freeing `rop->coeffs`.
void poly_clear(poly_t rop)
{
    if (rop == (struct _poly_t*) NULL)
        return;
    if (rop->coeffs == (mpfr_t*) NULL)
        return;
    
    // Clears each initialized coefficient
    for (size_t i = 0; i < rop->alloc; i++)
        mpfr_clear(rop->coeffs[i]);
    // Frees allocated array of coefficients
    free(rop->coeffs);

    // Sets coeffs to null pointer as a safeguard
    rop->coeffs = (mpfr_t*) NULL;
    rop->alloc = 0;
    rop->deg = 0;
}

// Variadic helper macro. Applies method to each argument.
#define DEFINE_POLY_AP(NAME, METHOD)                    \
    void poly_##NAME(poly_t first, ...)                 \
    {                                                   \
        va_list args;                                   \
        va_start(args, first);                          \
                                                        \
        struct _poly_t* current = first;                \
        while (current != (struct _poly_t*) NULL)       \
        {                                               \
            METHOD(current);                            \
            current = va_arg(args, struct _poly_t*);    \
        }                                               \
                                                        \
        va_end(args);                                   \
    }

DEFINE_POLY_AP(inits,   poly_init)
DEFINE_POLY_AP(clears,  poly_clear)

/* ************************************************************************* */

// Declares polynomial given array of coefficients and degree.
void poly_new_fr(poly_t rop, mpfr_t* coeffs, size_t deg)
{
    // TODO: document
    POLY_GUARD_INIT(ptr, scratch, rop,
        ptr_range(coeffs, rop->coeffs, rop->alloc));

    poly_set_deg(ptr, deg);
    for (size_t i = 0; i <= deg; i++)
        mpfr_set(ptr->coeffs[i], coeffs[i], MPFR_RNDN);
    poly_normalise(ptr);

    POLY_GUARD_CLEAR(rop, scratch);
}

// Declares polynomial given array of coefficients and degree.
#define DEFINE_POLY_NEW(SUFFIX, TYPE, SET_METHOD)                   \
    void poly_new_##SUFFIX(poly_t rop, TYPE* coeffs, size_t deg)    \
    {                                                               \
        poly_set_deg(rop, deg);                                     \
        for (size_t i = 0; i <= deg; i++)                           \
            SET_METHOD(rop->coeffs[i], coeffs[i], MPFR_RNDN);       \
        poly_normalise(rop);                                        \
    }

DEFINE_POLY_NEW(z,  mpz_t,                  mpfr_set_z)
DEFINE_POLY_NEW(ui, const unsigned long,    mpfr_set_ui)
DEFINE_POLY_NEW(si, const long,             mpfr_set_si)

/* ************************************************************************* */

// Sets polynomial `rop` to a copy of polynomial `op`.
void poly_set(poly_t rop, poly_t op)
{
    if (rop == op)
        return;

    // Sets degree of `op` to `rop`
    poly_set_deg(rop, op->deg);
    // Sets coefficients of `op` to `rop`
    for (size_t i = 0; i <= op->deg; i++)
        mpfr_set(rop->coeffs[i], op->coeffs[i], MPFR_RNDN);
}


// Sets `rop` to constant polynomial `rop(x) = op`, where `op` is a scalar.
void poly_set_fr(poly_t rop, mpfr_t op)
{
    // TODO: document
    POLY_SCALAR_GUARD_INIT(ptr, scratch, rop, op);

    poly_set_deg(rop, 0);
    mpfr_set(rop->coeffs[0], ptr, MPFR_RNDN);

    POLY_SCALAR_GUARD_CLEAR(scratch);
}

// Sets `rop` to constant polynomial `rop(x) = op`, where `op` is a scalar.
#define DEFINE_POLY_SET(SUFFIX, TYPE, SET_METHOD)       \
    void poly_set_##SUFFIX(poly_t rop, TYPE op)         \
    {                                                   \
        poly_set_deg(rop, 0);                           \
        SET_METHOD(rop->coeffs[0], op, MPFR_RNDN);      \
    }

DEFINE_POLY_SET(z,  mpz_t,          mpfr_set_z)
DEFINE_POLY_SET(ui, unsigned long,  mpfr_set_ui)
DEFINE_POLY_SET(si, long,           mpfr_set_si)

/* ************************************************************************* */