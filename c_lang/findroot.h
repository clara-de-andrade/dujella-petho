#ifndef FINDROOT_H
#define FINDROOT_H

#include <mpfr.h>
#include "polynom.h"

/* ************************************************************************* */

void findroot(mpfr_t rop, poly_t poly, mpfr_t init);

/* ************************************************************************* */

#endif
