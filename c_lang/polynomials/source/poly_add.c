#include "poly_ring.h"
#include "internal.h"

/* ************************************************************************* */

// Sets `rop` to the negation of `op`.
void poly_neg(poly_t rop, poly_t op)
{
    // deg(-op1) = max(op1)
    poly_set_deg(rop, op->deg);
    
    // Computes the negation of `op` coefficient by coefficient
    for (size_t i = 0; i <= op->deg; i++)
        mpfr_neg(rop->coeffs[i], op->coeffs[i], MPFR_RNDN);
}

// Sets `rop` to the sum of `op1` and `op2`.
void poly_add(poly_t rop, poly_t op1, poly_t op2)
{
    // deg(op1 + op2) <= max(deg(op1), deg(op2))
    size_t deg = max(op1->deg, op2->deg);

    poly_set_deg(rop, deg);
    // Computes the sum of `op1` and `op2` coefficient by coefficient
    for (size_t i = 0; i <= deg; i++)
        if (i <= op1->deg && i <= op2->deg)
            mpfr_add(rop->coeffs[i],
                op1->coeffs[i], op2->coeffs[i], MPFR_RNDN);
        else
        if (i <= op1->deg)
            mpfr_set(rop->coeffs[i], op1->coeffs[i], MPFR_RNDN);
        else
            mpfr_set(rop->coeffs[i], op2->coeffs[i], MPFR_RNDN);

    // Update `rop->deg` to true value
    poly_normalise(rop);
}

// Sets `rop` to the difference of `op1` and `op2`.
void poly_sub(poly_t rop, poly_t op1, poly_t op2)
{
    // deg(op1 - op2) <= max(deg(op1), deg(op2))
    size_t deg = max(op1->deg, op2->deg);
    
    poly_set_deg(rop, deg);
    // Computes the difference of `op1` and `op2` coefficient by coefficient
    for (size_t i = 0; i <= deg; i++)
        if (i <= op1->deg && i <= op2->deg)
            mpfr_sub(rop->coeffs[i],
                op1->coeffs[i], op2->coeffs[i], MPFR_RNDN);
        else
        if (i <= op1->deg)
            mpfr_set(rop->coeffs[i], op1->coeffs[i], MPFR_RNDN);
        else
            mpfr_neg(rop->coeffs[i], op2->coeffs[i], MPFR_RNDN);

    // Update `rop->deg` to true value
    poly_normalise(rop);
}

/* ************************************************************************* */