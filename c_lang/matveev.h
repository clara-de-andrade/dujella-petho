#ifndef MATVEEV_H
#define MATVEEV_H

#include <stddef.h>

#include <gmp.h>
#include <mpfr.h>

struct _mat_t {
    mpfr_prec_t prec;
    size_t m;

    size_t degree;
    mpfr_t* alpha;
    mpfr_t* A;

    mpfr_t C;
};

typedef struct _mat_t mat_t[1];

/*****************************************************************************/

void mat_init(mat_t rop, mpfr_prec_t prec, size_t m);
void mat_clear(mat_t rop);

void matveev_C(mat_t mat);

/*****************************************************************************/

#endif
