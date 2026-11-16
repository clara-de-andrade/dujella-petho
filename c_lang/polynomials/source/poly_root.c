#include "poly_ring.h"
#include "poly_root.h"
#include "internal.h"

/* ************************************************************************* */

// Sets `rop` to `eps * max(1, |x|)`.
static void set_tol(mpfr_t rop, mpfr_t op, mpfr_t eps)
{
    MPFR_GUARD_INIT(ptr, scratch, rop, rop == op);

    mpfr_abs(ptr, op, MPFR_RNDN);
    if (mpfr_cmp_ui(ptr, 1) < 0)
        mpfr_set_ui(ptr, 1, MPFR_RNDN);
    mpfr_mul(ptr, ptr, eps, MPFR_RNDN);

    MPFR_GUARD_CLEAR(scratch);
}

/* ************************************************************************* */

// Sets `rop` to root of `poly` closest to `init`. Returns 0 on convergence,
// positive if exceeded maximum number of iterations, and negative if either
// failed to convergence or found a stationary point. 
int poly_root(mpfr_t rop, poly_t poly, mpfr_t init)
{
    int status = POLY_ROOT_MAXITER;

    // TODO: document
    mpfr_prec_t prec = mpfr_get_prec(rop);
    size_t max_iter = 64 + 2 * (size_t) prec;

    mpfr_t x, x_new, P_eval, dP_eval, step, tol, eps;
    mpfr_inits2(prec,
        x, x_new, P_eval, dP_eval, step, tol, eps, (mpfr_ptr) NULL);
    
    mpfr_set(x, init, MPFR_RNDN);
    mpfr_set_ui_2exp(eps, 1, 1 - (mpfr_exp_t) prec, MPFR_RNDN);

    poly_t d_poly;
    poly_init(d_poly);
    poly_d(d_poly, poly);

    for (size_t iter = 0; iter < max_iter; iter++)
    {
        poly_eval_fr(P_eval, poly, x);
        if (mpfr_zero_p(P_eval))
        {
            status = POLY_ROOT_CONVERGENCE;
            break;
        }

        poly_eval_fr(dP_eval, d_poly, x);
        if (mpfr_zero_p(dP_eval))
        {
            status = POLY_ROOT_STATIONARY;
            break;
        }

        mpfr_div(step, P_eval, dP_eval, MPFR_RNDN);
        set_tol(tol, x, eps);
        if (mpfr_cmpabs(step, tol) <= 0)
        {
            status = POLY_ROOT_CONVERGENCE;
            break;
        }

        mpfr_sub(x_new, x, step, MPFR_RNDN);
        if (mpfr_cmp(x, x_new) == 0)
        {
            status = POLY_ROOT_CONVERGENCE;
            break;
        }
        mpfr_set(x, x_new, MPFR_RNDN);
    }

    mpfr_set(rop, x, MPFR_RNDN);
    mpfr_clears(x, x_new, P_eval, dP_eval, step, tol, eps, (mpfr_ptr) NULL);
    poly_clear(d_poly);

    return status;
}

/* ************************************************************************* */
