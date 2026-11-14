#include "findroot.h"

/* ************************************************************************* */

void findroot(mpfr_t rop, poly_t poly, mpfr_t init)
{
    mpfr_prec_t prec = mpfr_get_prec(rop);

    mpfr_t x, x_new, P_eval, dP_eval, step, tol, eps;
    mpfr_inits2(prec, x, x_new, P_eval, dP_eval, step, tol, eps, NULL);

    mpfr_init_set(x, init, MPFR_RNDN);
    mpfr_set_ui_2exp(eps, 1, 1 - (mpfr_exp_t) prec, MPFR_RNDN);

    poly_t d_poly;
    poly_init(d_poly);
    poly_d(d_poly, poly);

    while (1)
    {
        poly_eval_fr(P_eval, poly, x);
        if (mpfr_zero_p(P_eval))
            break;

        poly_eval_fr(dP_eval, d_poly, x);
        if (mpfr_zero_p(dP_eval))
            break;

        mpfr_div(step, P_eval, dP_eval, MPFR_RNDN);
        mpfr_abs(tol, x, MPFR_RNDN);
        mpfr_mul(tol, tol, eps, MPFR_RNDN);
        if (mpfr_cmpabs(step, tol) <= 0)
            break;

        mpfr_sub(x_new, x, step, MPFR_RNDN);
        if (mpfr_cmp(x, x_new) == 0)
            break;
        mpfr_set(x, x_new, MPFR_RNDN);
    }

    mpfr_set(rop, x, MPFR_RNDN);
    mpfr_clears(x, x_new, P_eval, dP_eval, step, tol, eps, NULL);
    poly_clear(d_poly);
}

/* ************************************************************************* */
