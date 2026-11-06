#ifndef __QKFIB_H__
#define __QKFIB_H__

#include <gmp.h>

/* ************************************************************************* */

// Computes the n-th term of the (q, k)-generalised Fibonacci sequence, as the
// n-th coefficient of the formal power series for its generating function,
// using the Bostan-Mori algorithm. (Algorithm 1 of arXiv:2008:08822v1.)
void qkfib(mpz_t op, mpz_t q, mpz_t k, mpz_t n);

/* ************************************************************************* */


#endif
