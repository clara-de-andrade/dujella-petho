#include "matveev.h"

#include <stdlib.h>

/*****************************************************************************/

void mat_init(mat_t rop, mpfr_prec_t prec, size_t m)
{
    rop->prec = prec;
    rop->m = m;

    rop->alpha = calloc(m, sizeof(mpfr_t));
    rop->A = calloc(m, sizeof(mpfr_t));
    for (size_t i = 0; i < m; i++)
        mpfr_inits2(prec, rop->alpha[i], rop->A[i], NULL);

    mpfr_init2(rop->C, prec);
}

void mat_clear(mat_t rop)
{
    for (size_t i = 0; i < rop->m; i++)
        mpfr_clears(rop->alpha[i], rop->A[i], NULL);
    free(rop->alpha);
    free(rop->A);

    mpfr_clear(rop->C);
}

/*****************************************************************************/

void matveev_C(mat_t mat)
{
    mpfr_t temp;
    mpfr_init2(temp, mat->prec);

    mpfr_set_str(mat->C, "1.4", 10, MPFR_RNDU);

    mpfr_set_ui(temp, 30, MPFR_RNDU);
    mpfr_pow_ui(temp, temp, mat->m + 3, MPFR_RNDU);
    mpfr_mul(mat->C, mat->C, temp, MPFR_RNDU);

    mpfr_set_ui(temp, 3, MPFR_RNDU);
    mpfr_pow_ui(temp, temp, 4, MPFR_RNDU);
    mpfr_mul(mat->C, mat->C, temp, MPFR_RNDU);
    mpfr_sqrt_ui(temp, mat->m, MPFR_RNDU);
    mpfr_mul(mat->C, mat->C, temp, MPFR_RNDU);

    mpfr_set_ui(temp, mat->degree, MPFR_RNDU);
    mpfr_pow_ui(temp, temp, 2, MPFR_RNDU);
    mpfr_mul(mat->C, mat->C, temp, MPFR_RNDU);

    mpfr_set_ui(temp, mat->degree, MPFR_RNDU);
    mpfr_log(temp, temp, MPFR_RNDU);
    mpfr_add_ui(temp, temp, 1, MPFR_RNDU);
    mpfr_mul(mat->C, mat->C, temp, MPFR_RNDU);

    for (size_t i = 0; i < mat->m; i++)
        mpfr_mul(mat->C, mat->C, mat->A[i], MPFR_RNDU);

    mpfr_clear(temp);
}

/*****************************************************************************/
