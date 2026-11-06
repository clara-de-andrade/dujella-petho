#ifndef __RECSEQ_H__
#define __RECSEQ_H__

#include <gmp.h>

/* ************************************************************************* */

// Computes the n-th term of a linear recurrent sequence, as the n-th
// coefficient of the formal power series for its generating function,
// using the Bostan-Mori algorithm. (Algorithm 1 of arXiv:2008.08822v1.)
void reqseq(mpz_t op, mpz_t k, mpz_t* coeffs, mpz_t* inits, mpz_t n);

/* ************************************************************************* */


#endif
