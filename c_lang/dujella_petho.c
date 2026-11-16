#include "dujella_petho.h"

/* ************************************************************************* */

// Initializes parameters of `rop` with precision `prec`.
void dp_init(dp_t rop, mpfr_prec_t prec)
{
    rop->prec = prec;
    mpfr_inits2(prec, rop->gamma, rop->mu, rop->A, rop->B, rop->eps, NULL);
    mpz_inits(rop->M, rop->q, rop->N, NULL);
}

// Clears parameters of `rop`.
void dp_clear(dp_t rop)
{
    mpfr_clears(rop->gamma, rop->mu, rop->A, rop->B, rop->eps, NULL);
    mpz_clears(rop->M, rop->q, rop->N, NULL);
}

/* ************************************************************************* */

// Computes dp->q, dp->eps and dp->N from remaining parameters.
void dujella_petho(dp_t dp)
{
    mpz_t target;
    mpz_init(target);
    mpz_mul_ui(target, dp->M, 6);

    mpz_t q, q_old, q_new, a;
    mpz_inits(q, q_old, q_new, a, NULL);

    mpz_set_ui(q, 0);
    mpz_set_ui(q_old, 1);

    mpfr_t eta, temp, t1, t2, eps;
    mpfr_inits2(dp->prec, eta, temp, t1, t2, eps, NULL);

    mpfr_set(eta, dp->gamma, MPFR_RNDN);

    while (1)
    {
        mpfr_get_z(a, eta, MPFR_RNDD);
        mpfr_sub_z(eta, eta, a, MPFR_RNDN);

        mpz_mul(q_new, q, a);
        mpz_add(q_new, q_new, q_old);
        mpz_set(q_old, q);
        mpz_set(q, q_new);

        if (mpfr_zero_p(eta))
            break;
        mpfr_ui_div(eta, 1, eta, MPFR_RNDN);

        if (mpz_cmp(q, target) <= 0)
            continue;

        mpfr_mul_z(t1, dp->mu, q, MPFR_RNDN);
        mpfr_round(temp, t1);
        mpfr_sub(t1, t1, temp, MPFR_RNDN);
        mpfr_abs(t1, t1, MPFR_RNDN);

        mpfr_mul_z(t2, dp->gamma, q, MPFR_RNDN);
        mpfr_round(temp, t2);
        mpfr_sub(t2, t2, temp, MPFR_RNDN);
        mpfr_abs(t2, t2, MPFR_RNDN);
        mpfr_mul_z(t2, t2, dp->M, MPFR_RNDN);

        mpfr_sub(eps, t1, t2, MPFR_RNDN);
        if (mpfr_sgn(eps) <= 0)
            continue;

        break;
    }

    mpz_set(dp->q, q);
    mpfr_set(dp->eps, eps, MPFR_RNDN);

    mpz_clears(target, q, q_old, q_new, a, NULL);
    mpfr_clears(eta, t1, t2, eps, NULL);

    mpfr_t num, den;
    mpfr_inits2(dp->prec, num, den, NULL);

    mpfr_mul_z(num, dp->A, dp->q, MPFR_RNDN);
    mpfr_div(num, num, dp->eps, MPFR_RNDN);
    mpfr_log(num, num, MPFR_RNDN);
    mpfr_log(den, dp->B, MPFR_RNDN);
    mpfr_div(temp, num, den, MPFR_RNDN);
    
    mpfr_get_z(dp->N, temp, MPFR_RNDD);

    mpfr_clears(num, den, temp, NULL);
}

/* ************************************************************************* */