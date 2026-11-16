#include "internal.h"

/* ************************************************************************* */

// Reallocates `rop->coeffs` to either twice as many or `alloc` elements.
void poly_realloc(poly_t rop, size_t alloc)
{
    // Does nothing if `rop->coeffs` has enough elements allocated
    if (alloc <= rop->alloc)
        return;

    size_t max_alloc = SIZE_MAX / sizeof(mpfr_t);

    // OVERFLOW ERROR
    if (alloc > max_alloc)
        exit(EXIT_FAILURE);

    size_t new_alloc;
    if (rop->alloc > (max_alloc >> 1))
        // If `rop->alloc` is too large, allocate the max number of elements
        new_alloc = max_alloc;
    else
        // Otherwise, allocate either twice as many or `alloc` elements
        new_alloc = max(alloc, rop->alloc << 1);
    mpfr_t* coeffs = realloc(rop->coeffs, new_alloc * sizeof(mpfr_t));

    // MEMORY ERROR
    if (coeffs == (mpfr_t*) NULL)
        exit(EXIT_FAILURE);
    rop->coeffs = coeffs;

    for (size_t i = rop->alloc; i < new_alloc; i++)
        // Initializes coefficients and sets them to zero by default
        mpfr_init_set_ui(rop->coeffs[i], 0, MPFR_RNDN);
    rop->alloc = new_alloc;
}

// Sets `rop->deg` to `deg`. If needed, reallocates or truncates `rop->coeff`.
void poly_set_deg(poly_t rop, size_t deg)
{
    // OVERFLOW ERROR
    if (deg == SIZE_MAX)
        exit(EXIT_FAILURE);

    // Reallocates `rop` to at least `deg + 1` elements
    poly_realloc(rop, deg + 1);

    // If `deg < rop->deg`, set coefficients of higher degree to zero
    if (deg < rop->deg)
        for (size_t i = deg + 1; i <= rop->deg; i++)
            mpfr_set_ui(rop->coeffs[i], 0, MPFR_RNDN);

    // Set `rop->deg` to `deg` (not necessarily the true degree!)
    rop->deg = deg;
}

/* ************************************************************************* */

// Sets `rop->deg` to largest `index` such that `rop->coeff[index] != 0`.
void poly_normalise(poly_t rop)
{
    while (rop->deg > 0)
    {
        if (mpfr_zero_p(rop->coeffs[rop->deg]))
            rop->deg--;
        else
            break;
    }
}

/* ************************************************************************* */