#include "polynom.h"

#include <stdio.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdlib.h>

#include <string.h>

/* ************************************************************************* */

#define POLY_MUL_CUTOFF 32

#define do_nothing(x) (x)
#define max(a, b) ((a) >= (b) ? (a) : (b))
#define min(a, b) ((a) <= (b) ? (a) : (b))

/* ************************************************************************* */

void _poly_ap(void (*method)(struct _poly_t*), struct _poly_t* first, ...)
{
    va_list args;
    va_start(args, first);

    struct _poly_t* current = first;
    while (current != (struct _poly_t*) NULL)
    {
        method(current);
        current = va_arg(args, struct _poly_t*);
    }
    
    va_end(args);
}

void poly_init(poly_t rop)
{
    rop->coeffs = malloc(sizeof(mpfr_t));
    if (rop->coeffs == (mpfr_t*) NULL)
        exit(EXIT_FAILURE);
    mpfr_init_set_ui(rop->coeffs[0], 0, MPFR_RNDN);

    rop->alloc = 1;
    rop->deg = 0;
}

void poly_clear(poly_t rop)
{
    if (rop == (struct _poly_t*) NULL)
        return;

    if (rop->coeffs == (mpfr_t*) NULL)
        return;
    
    for (size_t i = 0; i < rop->alloc; i++)
        mpfr_clear(rop->coeffs[i]);
    free(rop->coeffs);

    rop->coeffs = (mpfr_t*) NULL;
}

/* ************************************************************************* */

static void poly_realloc(poly_t rop, size_t alloc)
{
    if (alloc <= rop->alloc)
        return;

    if (alloc > (SIZE_MAX / sizeof(mpfr_t)))
        exit(EXIT_FAILURE);

    size_t new_alloc;
    if (alloc > (SIZE_MAX / sizeof(mpfr_t)) / 2)
        new_alloc = alloc;
    else
        new_alloc = alloc << 1;

    mpfr_t* coeffs = realloc(rop->coeffs, new_alloc * sizeof(mpfr_t));
    if (coeffs == (mpfr_t*) NULL)
        exit(EXIT_FAILURE);
    rop->coeffs = coeffs;

    for (size_t i = rop->alloc; i < new_alloc; i++)
        mpfr_init_set_ui(rop->coeffs[i], 0, MPFR_RNDN);
    rop->alloc = new_alloc;
}

void poly_set(poly_t rop, poly_t op)
{
    poly_set_deg(rop, op->deg);

    for (size_t i = 0; i <= op->deg; i++)
        mpfr_set(rop->coeffs[i], op->coeffs[i], MPFR_RNDN);
}

#define DEFINE_POLY_SET(SUFFIX, TYPE, SET_METHOD)   \
    void poly_set_##SUFFIX(poly_t rop, TYPE op)     \
    {                                               \
        poly_set_deg(rop, 0);                       \
        SET_METHOD(rop->coeffs[0], op, MPFR_RNDN);  \
    }

DEFINE_POLY_SET(fr, mpfr_t,         mpfr_set)
DEFINE_POLY_SET(z,  mpz_t,          mpfr_set_z)
DEFINE_POLY_SET(ui, unsigned long,  mpfr_set_ui)
DEFINE_POLY_SET(si, long,           mpfr_set_si)

/* ************************************************************************* */

static void poly_normalise(poly_t rop)
{
    while (rop->deg > 0)
    {
        if (mpfr_zero_p(rop->coeffs[rop->deg]))
            rop->deg--;
        else
            break;
    }
}

size_t poly_get_deg(poly_t op)
{
    return op->deg;
}

void poly_set_deg(poly_t rop, size_t deg)
{
    if (deg == SIZE_MAX)
        exit(EXIT_FAILURE);
    poly_realloc(rop, deg + 1);

    if (deg < rop->deg)
        for (size_t i = deg + 1; i <= rop->deg; i++)
            mpfr_set_ui(rop->coeffs[i], 0, MPFR_RNDN);
    else
    if (deg > rop->deg)
        for (size_t i = rop->deg + 1; i <= deg; i++)
            mpfr_set_ui(rop->coeffs[i], 0, MPFR_RNDN);

    rop->deg = deg;
}


void poly_get_coeff(mpfr_t rop, poly_t op, size_t index)
{
    if (index > op->deg)
        mpfr_set_ui(rop, 0, MPFR_RNDN);
    else
        mpfr_set(rop, op->coeffs[index], MPFR_RNDN);
}

#define DEFINE_POLY_SET_COEFF(SUFFIX, TYPE, SET_METHOD, IS_ZERO)        \
    void poly_set_coeff_##SUFFIX(poly_t rop, TYPE op, size_t index)     \
    {                                                                   \
        int op_is_zero = (IS_ZERO(op) == 0);                            \
                                                                        \
        if (index > rop->deg)                                           \
        {                                                               \
            if (op_is_zero)                                             \
                return;                                                 \
                                                                        \
            poly_set_deg(rop, index);                                   \
            SET_METHOD(rop->coeffs[index], op, MPFR_RNDN);              \
        }                                                               \
        else                                                            \
        {                                                               \
            SET_METHOD(rop->coeffs[index], op, MPFR_RNDN);              \
            if (index == rop->deg && op_is_zero)                        \
                poly_normalise(rop);                                    \
        }                                                               \
    }

DEFINE_POLY_SET_COEFF(fr, mpfr_t,         mpfr_set,     mpfr_sgn)
DEFINE_POLY_SET_COEFF(z,  mpz_t,          mpfr_set_z,   mpz_sgn)
DEFINE_POLY_SET_COEFF(ui, unsigned long,  mpfr_set_ui,  do_nothing)
DEFINE_POLY_SET_COEFF(si, long,           mpfr_set_si,  do_nothing)

/* ************************************************************************* */

void poly_eval(poly_t rop, poly_t op1, poly_t op2)
{
    if (op1->deg == 0)
    {
        mpfr_t coeff;
        mpfr_init_set(coeff, op1->coeffs[0], MPFR_RNDN);

        poly_set_fr(rop, coeff);
        poly_normalise(rop);

        mpfr_clear(coeff);
        return;
    }

    if (op2->deg == 0)
    {
        mpfr_t coeff;
        mpfr_init(coeff);

        poly_eval_fr(coeff, op1, op2->coeffs[0]);
        poly_set_fr(rop, coeff);

        mpfr_clear(coeff);
        return;
    }

    size_t d = op1->deg + 1;
    if (d > (SIZE_MAX >> 1) + 1)
        exit(EXIT_FAILURE);
    d = 1ULL << (sizeof(unsigned long long) * 8 - __builtin_clzll(d - 1));

    size_t m = d >> 1;

    poly_t P0, P1;
    poly_inits(P0, P1);

    for (size_t i = 0; i <= op1->deg; i++)
        if (i < m)
            poly_set_coeff_fr(P0, op1->coeffs[i], i);
        else
            poly_set_coeff_fr(P1, op1->coeffs[i], i - m);

    poly_t P0_eval, P1_eval;
    poly_inits(P0_eval, P1_eval);

    poly_eval(P0_eval, P0, op2);
    poly_eval(P1_eval, P1, op2);

    poly_t temp;
    poly_init(temp);

    poly_pow_ui(temp, op2, m);
    poly_mul(P1_eval, P1_eval, temp);
    poly_add(rop, P0_eval, P1_eval);

    poly_clears(P0, P0_eval, P1, P1_eval, temp);
}

#define DEFINE_POLY_EVAL(SUFFIX, TYPE, INIT_SET_METHOD)                         \
    void poly_eval_##SUFFIX(mpfr_t rop, poly_t op1, TYPE op2)                   \
    {                                                                           \
        mpfr_t acc;                                                             \
        mpfr_init_set_ui(acc, 0, MPFR_RNDN);                                    \
                                                                                \
        mpfr_t op2_fr;                                                          \
        INIT_SET_METHOD(op2_fr, op2, MPFR_RNDN);                                \
                                                                                \
        for (size_t i = 0; i <= op1->deg; i++)                                  \
            mpfr_fma(acc, acc, op2_fr, op1->coeffs[op1->deg - i], MPFR_RNDN);   \
        mpfr_set(rop, acc, MPFR_RNDN);                                          \
                                                                                \
        mpfr_clears(acc, op2_fr, NULL);                                         \
    }

DEFINE_POLY_EVAL(fr,    mpfr_t,         mpfr_init_set)
DEFINE_POLY_EVAL(z,     mpz_t,          mpfr_init_set_z)
DEFINE_POLY_EVAL(ui,    unsigned long,  mpfr_init_set_ui)
DEFINE_POLY_EVAL(si,    long,           mpfr_init_set_si)

/* ************************************************************************* */

static void str_append(char** buffer, size_t* len, size_t* cap, char* str)
{
    size_t add = strlen(str);
    if (*len + add + 1 > *cap)
    {
        while (*len + add + 1 > *cap)
            *cap <<= 1;
        char* temp = realloc(*buffer, *cap);
        if (temp == NULL)
            exit(EXIT_FAILURE);
        *buffer = temp;
    }
    memcpy(*buffer + *len, str, add + 1);
    *len += add;
        
}

char* poly_str(poly_t op)
{
    char* coeff_str = NULL;

    if (op->deg == 0)
    {
        if (mpfr_asprintf(&coeff_str, "%Rg", op->coeffs[0]) < 0)
            exit(EXIT_FAILURE);
        return coeff_str;
    }

    size_t cap = 64;
    size_t len = 0;
    char* res = malloc(cap);
    if (res == NULL)
        exit(EXIT_FAILURE);
    res[0] = '\0';

    mpfr_t abs_coeff;
    mpfr_init(abs_coeff);

    int has_terms = 0;
    char exp_buf[32];

    for (size_t i = op->deg; i != (size_t) -1; i--)
    {
        if (mpfr_zero_p(op->coeffs[i]))
            continue;

        int is_neg = (mpfr_sgn(op->coeffs[i]) < 0);
        mpfr_abs(abs_coeff, op->coeffs[i], MPFR_RNDN);

        if (has_terms)
            str_append(&res, &len, &cap, is_neg ? " - " : " + ");
        else
        if (is_neg)
            str_append(&res, &len, &cap, "-");

        if (!(i > 0 && mpfr_cmp_ui(abs_coeff, 1) == 0))
        {
            if (mpfr_asprintf(&coeff_str, "%Rg", abs_coeff) < 0)
                exit(EXIT_FAILURE);
            str_append(&res, &len, &cap, coeff_str);
            mpfr_free_str(coeff_str);
        }

        if (i == 1)
            str_append(&res, &len, &cap, "x");
        else
        if (i > 1)
        {
            snprintf(exp_buf, sizeof(exp_buf), "x**%zu", i);
            str_append(&res, &len, &cap, exp_buf);
        }

        has_terms = 1;
    }

    mpfr_clear(abs_coeff);

    if (!has_terms)
        str_append(&res, &len, &cap, "0");

    return res;
}

/* ************************************************************************* */

void poly_neg(poly_t rop, poly_t op)
{
    poly_set_deg(rop, op->deg);

    for (size_t i = 0; i <= op->deg; i++)
        mpfr_neg(rop->coeffs[i], op->coeffs[i], MPFR_RNDN);

    poly_normalise(rop);
}

void poly_add(poly_t rop, poly_t op1, poly_t op2)
{
    size_t deg = max(op1->deg, op2->deg);
    
    poly_set_deg(rop, deg);

    for (size_t i = 0; i <= deg; i++)
        if (i <= op1->deg && i <= op2->deg)
            mpfr_add(rop->coeffs[i], op1->coeffs[i], op2->coeffs[i], MPFR_RNDN);
        else if (i <= op1->deg)
            mpfr_set(rop->coeffs[i], op1->coeffs[i], MPFR_RNDN);
        else
            mpfr_set(rop->coeffs[i], op2->coeffs[i], MPFR_RNDN);

    poly_normalise(rop);
}

void poly_sub(poly_t rop, poly_t op1, poly_t op2)
{
    size_t deg = max(op1->deg, op2->deg);
    
    poly_set_deg(rop, deg);

    for (size_t i = 0; i <= deg; i++)
        if (i <= op1->deg && i <= op2->deg)
            mpfr_sub(rop->coeffs[i], op1->coeffs[i], op2->coeffs[i], MPFR_RNDN);
        else if (i <= op1->deg)
            mpfr_set(rop->coeffs[i], op1->coeffs[i], MPFR_RNDN);
        else
            mpfr_neg(rop->coeffs[i], op2->coeffs[i], MPFR_RNDN);

    poly_normalise(rop);
}

/* ************************************************************************* */

#define poly_iszero(op) (op->deg == 0 && mpfr_zero_p(op->coeffs[0]))

static void poly_cauchy(poly_t rop, poly_t op1, poly_t op2)
{
    size_t deg = op1->deg + op2->deg;

    poly_t res;
    poly_init(res);
    poly_set_deg(res, deg);

    mpfr_t cauchy;
    mpfr_init(cauchy);

    for (size_t k = 0; k <= deg; k++)
    {
        mpfr_set_ui(cauchy, 0, MPFR_RNDN);

        size_t i_min = (k > op2->deg) ? (k - op2->deg) : 0;
        size_t i_max = min(k, op1->deg);

        for (size_t i = i_min; i <= i_max; i++)
        {
            if (mpfr_sgn(op1->coeffs[i]) == 0 || mpfr_sgn(op2->coeffs[k - i]) == 0)
                continue;
            mpfr_fma(cauchy, op1->coeffs[i], op2->coeffs[k - i], cauchy, MPFR_RNDN);
        }
        mpfr_set(res->coeffs[k], cauchy, MPFR_RNDN);
    }

    poly_normalise(res);
    poly_set(rop, res);

    poly_clear(res);
    mpfr_clear(cauchy);
}

void poly_mul(poly_t rop, poly_t op1, poly_t op2)
{
    if (poly_iszero(op1) || poly_iszero(op2))
    {
        poly_set_ui(rop, 0);
        return;
    }

    if (op1->deg > SIZE_MAX - op2->deg)
        exit(EXIT_FAILURE);

    size_t d = max(op1->deg, op2->deg) + 1;
    if (d > (SIZE_MAX >> 1) + 1)
        exit(EXIT_FAILURE);

    if (d == 1)
    {
        mpfr_t coeff;
        mpfr_init(coeff);

        mpfr_mul(coeff, op1->coeffs[0], op2->coeffs[0], MPFR_RNDN);
        poly_set_fr(rop, coeff);

        mpfr_clear(coeff);
        return;
    }

    if (d <= POLY_MUL_CUTOFF)
    {
        poly_cauchy(rop, op1, op2);
        return;
    }
    d = 1ULL << (sizeof(unsigned long long) * 8 - __builtin_clzll(d - 1));

    size_t m = d >> 1;

    poly_t P0, P1, Q0, Q1;
    poly_inits(P0, P1, Q0, Q1);

    for (size_t i = 0; i <= op1->deg; i++)
        if (i < m)
            poly_set_coeff_fr(P0, op1->coeffs[i], i);
        else
            poly_set_coeff_fr(P1, op1->coeffs[i], i - m);

    for (size_t i = 0; i <= op2->deg; i++)
        if (i < m)
            poly_set_coeff_fr(Q0, op2->coeffs[i], i);
        else
            poly_set_coeff_fr(Q1, op2->coeffs[i], i - m);

    poly_t z0, z1, z2, temp, middle;
    poly_inits(z0, z1, z2, temp, middle);

    poly_mul(z0, P0, Q0);
    poly_mul(z2, P1, Q1);

    poly_add(z1, P0, P1);
    poly_add(temp, Q0, Q1);
    poly_mul(z1, z1, temp);

    poly_add(middle, z0, z2);
    poly_sub(middle, z1, middle);

    poly_set(rop, z0);
    poly_shift(temp, middle, m);
    poly_add(rop, rop, temp);
    poly_shift(temp, z2, d);
    poly_add(rop, rop, temp);

    poly_clears(P0, P1, Q0, Q1);
    poly_clears(z0, z1, z2, temp, middle);
}

#define DEFINE_POLY_MUL(SUFFIX, TYPE, MUL_METHOD, IS_ZERO)                  \
    void poly_mul_##SUFFIX(poly_t rop, poly_t op1, TYPE op2)                \
    {                                                                       \
        if (poly_iszero(op1) || (IS_ZERO(op2) == 0))                        \
        {                                                                   \
            poly_set_ui(rop, 0);                                            \
            return;                                                         \
        }                                                                   \
                                                                            \
        poly_set_deg(rop, op1->deg);                                        \
                                                                            \
        for (size_t i = 0; i <= op1->deg; i++)                              \
            MUL_METHOD(rop->coeffs[i], op1->coeffs[i], op2, MPFR_RNDN);     \
                                                                            \
        poly_normalise(rop);                                                \
    }

DEFINE_POLY_MUL(fr,     mpfr_t,         mpfr_mul,       mpfr_sgn)
DEFINE_POLY_MUL(z,      mpz_t,          mpfr_mul_z,     mpz_sgn)
DEFINE_POLY_MUL(ui,     unsigned long,  mpfr_mul_ui,    do_nothing)
DEFINE_POLY_MUL(si,     long,           mpfr_mul_si,    do_nothing)

void poly_shift(poly_t rop, poly_t op, size_t exp)
{
    if (poly_iszero(op))
    {
        poly_set_ui(rop, 0);
        return;
    }

    size_t old_deg = op->deg;
    if (exp > SIZE_MAX - old_deg)
        exit(EXIT_FAILURE);
    poly_set_deg(rop, old_deg + exp);

    for (size_t i = old_deg; i != (size_t) -1; i--)
        mpfr_set(rop->coeffs[i + exp], op->coeffs[i], MPFR_RNDN);

    for (size_t i = 0; i < exp; i++)
        mpfr_set_ui(rop->coeffs[i], 0, MPFR_RNDN);

    poly_normalise(rop);
}

/* ************************************************************************* */

void poly_pow_z(poly_t rop, poly_t op1, mpz_t op2)
{
    if (mpz_cmp_ui(op2, 0) == 0)
    {
        poly_set_ui(rop, 1);
        return;
    }

    if (mpz_cmp_ui(op2, 1) == 0)
    {
        poly_set(rop, op1);
        return;
    }

    poly_t acc, base;
    poly_init_set_ui(acc, 1);
    poly_init_set(base, op1);

    mpz_t exp;
    mpz_init_set(exp, op2);

    while (mpz_sgn(exp) > 0)
    {
        if (mpz_odd_p(exp))
            poly_mul(acc, acc, base);
        mpz_div_2exp(exp, exp, 1);
        if (mpz_cmp_ui(exp, 0) > 0)
            poly_mul(base, base, base);
    }

    poly_set(rop, acc);

    mpz_clear(exp);
    poly_clears(acc, base);
}

void poly_pow_ui(poly_t rop, poly_t op1, unsigned long op2)
{
    if (op2 == 0)
    {
        poly_set_ui(rop, 1);
        return;
    }

    if (op2 == 1)
    {
        poly_set(rop, op1);
        return;
    }

    poly_t acc, base;
    poly_init_set_ui(acc, 1);
    poly_init_set(base, op1);

    while (op2 > 0)
    {
        if (op2 & 1)
            poly_mul(acc, acc, base);
        op2 >>= 1;
        if (op2 > 0)
            poly_mul(base, base, base);
    }

    poly_set(rop, acc);

    poly_clears(acc, base);
}

/* ************************************************************************* */

void poly_d(poly_t rop, poly_t op)
{
    if (op->deg == 0)
    {
        poly_set_ui(rop, 0);
        return;
    }

    size_t orig_deg = op->deg;
    size_t old_deg = rop->deg;
    poly_realloc(rop, orig_deg);

    for (size_t i = 1; i <= orig_deg; i++)
        mpfr_mul_ui(rop->coeffs[i - 1], op->coeffs[i], i, MPFR_RNDN);

    for (size_t i = orig_deg; i <= old_deg; i++)
        mpfr_set_ui(rop->coeffs[i], 0, MPFR_RNDN);

    rop->deg = orig_deg - 1;
    poly_normalise(rop);
}

/* ************************************************************************* */