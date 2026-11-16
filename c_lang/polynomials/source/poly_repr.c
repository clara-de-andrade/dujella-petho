#include "internal.h"

#include <stdio.h>
#include <string.h>

/* ************************************************************************* */

// Static helper function. Appends `str`  to `buffer` of length `len` with
// `cap` allocated `char`s.
static void str_append(char** buffer, size_t* len, size_t* cap, const char* str)
{
    size_t add = strlen(str);
    if (*len + add + 1 > *cap)
    {
        while (*len + add + 1 > *cap)
            *cap <<= 1;
        char* temp = realloc(*buffer, *cap);
        // MEMORY ERROR
        if (temp == (char*) NULL)
            exit(EXIT_FAILURE);
        *buffer = temp;
    }
    memcpy(*buffer + *len, str, add + 1);
    *len += add;    
}

// Returns string representation of polynomial `op`.
char* poly_repr(poly_t op)
{
    // TODO: document
    char* coeff_str = (char*) NULL;
    size_t cap = 64;
    size_t len = 0;

    char* res = malloc(cap);
    // MEMORY ERROR
    if (res == (char*) NULL)
        exit(EXIT_FAILURE);
    res[0] = '\0';

    if (op->deg == 0)
    {
        // MEMORY ERROR
        if (mpfr_asprintf(&coeff_str, "%Rg", op->coeffs[0]) < 0)
            exit(EXIT_FAILURE);
        str_append(&res, &len, &cap, coeff_str);
        mpfr_free_str(coeff_str);
        return res;
    }

    mpfr_t abs_coeff;
    mpfr_init(abs_coeff);

    int has_terms = 0;
    char exp_buf[32];

    for (size_t i = op->deg; i != (size_t) -1; i--)
    {
        if (mpfr_zero(op->coeffs[i]))
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
            // MEMORY ERROR
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
