#ifndef DPTEST_H
#define DPTEST_H

#include <gmp.h>
#include <mpfr.h>

/* ************************************************************************* */

struct _dp_t {
    mpfr_prec_t prec;
    
    mpfr_t gamma;
    mpfr_t mu;
    mpfr_t A;
    mpfr_t B;
    mpfr_t epsilon;

    mpz_t M;
    mpz_t q;
    mpz_t N;
};
typedef struct _dp_t dp_t[1];

/* ************************************************************************* */

void dp_init(dp_t rop, mpfr_prec_t prec);
void dp_clear(dp_t rop);

/* ************************************************************************* */

void dptest(dp_t dp);

/* ************************************************************************* */

#endif
