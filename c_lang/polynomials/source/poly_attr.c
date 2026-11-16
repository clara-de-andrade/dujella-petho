#include "poly_type.h"
#include "internal.h"

/* ************************************************************************* */

// Returns `rop->deg`.
size_t poly_get_deg(poly_t op)
{
    return op->deg;
}

// Sets `rop` to `op->coeffs[index]`, setting to `0` if `index > op->deg`.
void poly_get_coeff(mpfr_t rop, poly_t op, size_t index)
{
    // If `index > op->deg`, set `rop` to zero
    if (index > op->deg)
        mpfr_set_ui(rop, 0, MPFR_RNDN);
    // Otherwise, set `rop` to `op->coeffs[index]`
    else
        mpfr_set(rop, op->coeffs[index], MPFR_RNDN);
}

/* ************************************************************************* */

// Sets `rop->coeffs[index]` to `op`, reallocating `rop` if `index > op->deg`.
void poly_set_coeff_fr(poly_t rop, mpfr_t op, size_t index)
{
    // TODO: document
    int op_is_zero = mpfr_zero_p(op);

    if (index > rop->deg)
    {
        if (op_is_zero)
            return;

        POLY_SCALAR_GUARD_INIT(ptr, scratch, rop, op);

        poly_set_deg(rop, index);
        mpfr_set(rop->coeffs[index], ptr, MPFR_RNDN);

        POLY_SCALAR_GUARD_CLEAR(scratch);
    }
    else
    {

        mpfr_set(rop->coeffs[index], op, MPFR_RNDN);
        if (index == rop->deg && op_is_zero)
            poly_normalise(rop);
    }
}

// Sets `rop->coeffs[index]` to `op`, reallocating `rop` if `index > op->deg`.
#define DEFINE_POLY_SET_COEFF(SUFFIX, TYPE, SET_METHOD, IS_ZERO)        \
    void poly_set_coeff_##SUFFIX(poly_t rop, TYPE op, size_t index)     \
    {                                                                   \
        int op_is_zero = IS_ZERO(op);                                   \
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

DEFINE_POLY_SET_COEFF(z,  mpz_t,          mpfr_set_z,   mpz_iszero)
DEFINE_POLY_SET_COEFF(ui, unsigned long,  mpfr_set_ui,  ui_iszero)
DEFINE_POLY_SET_COEFF(si, long,           mpfr_set_si,  si_iszero)

/* ************************************************************************* */