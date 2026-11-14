#ifndef POLYNOM_H
#define POLYNOM_H

#include <gmp.h>
#include <mpfr.h>

/* ************************************************************************* */

struct _poly_t {
    mpfr_t* coeffs;
    size_t alloc;
    size_t deg;
};

typedef struct _poly_t poly_t[1];

/* ************************************************************************* */

void _poly_ap(void (*method)(struct _poly_t*), struct _poly_t* first, ...);

void poly_init(poly_t rop);
#define poly_inits(...)  \
    _poly_ap(poly_init, __VA_ARGS__, (struct _poly_t*) NULL)

void poly_clear(poly_t rop);
#define poly_clears(...) \
    _poly_ap(poly_clear, __VA_ARGS__, (struct _poly_t*) NULL)

/* ************************************************************************* */

void poly_set(poly_t rop, poly_t op);
#define poly_init_set(rop, op)      \
    do {                            \
        poly_init(rop);             \
        poly_set(rop, op);          \
    } while (0)

void poly_set_fr(poly_t rop, mpfr_t op);
#define poly_init_set_fr(rop, op)   \
    do {                            \
        poly_init(rop);             \
        poly_set_fr(rop, op);       \
    } while (0)

void poly_set_z(poly_t rop, mpz_t op);
#define poly_init_set_z(rop, op)    \
    do {                            \
        poly_init(rop);             \
        poly_set_z(rop, op);        \
    } while (0)

void poly_set_ui(poly_t rop, unsigned long op);
#define poly_init_set_ui(rop, op)   \
    do {                            \
        poly_init(rop);             \
        poly_set_ui(rop, op);       \
    } while (0)

void poly_set_si(poly_t rop, long op);
#define poly_init_set_si(rop, op)   \
    do {                            \
        poly_init(rop);             \
        poly_set_si(rop, op);       \
    } while (0)

/* ************************************************************************* */

size_t poly_get_deg(poly_t rop);
void poly_set_deg(poly_t rop, size_t deg);

void poly_get_coeff(mpfr_t rop, poly_t op, size_t index);
void poly_set_coeff_fr(poly_t rop, mpfr_t op, size_t index);

void poly_set_coeff_z(poly_t rop, mpz_t op, size_t index);
void poly_set_coeff_ui(poly_t rop, unsigned long op, size_t index);
void poly_set_coeff_si(poly_t rop, long op, size_t index);

/* ************************************************************************* */

void poly_eval(poly_t rop, poly_t op1, poly_t op2);
void poly_eval_fr(mpfr_t rop, poly_t op1, mpfr_t op2);

void poly_eval_z(mpfr_t rop, poly_t op1, mpz_t op2);
void poly_eval_ui(mpfr_t rop, poly_t op1, unsigned long op2);
void poly_eval_si(mpfr_t rop, poly_t op1, long op2);

/* ************************************************************************* */

char* poly_str(poly_t op);

/* ************************************************************************* */

void poly_neg(poly_t rop, poly_t op);
void poly_add(poly_t rop, poly_t op1, poly_t op2);
void poly_sub(poly_t rop, poly_t op1, poly_t op2);

/* ************************************************************************* */

void poly_mul(poly_t rop, poly_t op1, poly_t op2);
void poly_mul_fr(poly_t rop, poly_t op1, mpfr_t op2);

void poly_mul_z(poly_t rop, poly_t op1, mpz_t op2);
void poly_mul_ui(poly_t rop, poly_t op1, unsigned long op2);
void poly_mul_si(poly_t rop, poly_t op1, long op2);

void poly_shift(poly_t rop, poly_t op, size_t exp);

/* ************************************************************************* */

void poly_pow_z(poly_t rop, poly_t op1, mpz_t op2);
void poly_pow_ui(poly_t rop, poly_t op1, unsigned long op2);

/* ************************************************************************* */

void poly_d(poly_t rop, poly_t op);

/* ************************************************************************* */

#endif
