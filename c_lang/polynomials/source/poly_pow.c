#include "poly_ring.h"
#include "internal.h"

/* ************************************************************************* */

// Sets `rop` to the power of `op1` by a non-negative integer `op2`.
void poly_pow_z(poly_t rop, poly_t op1, mpz_t op2)
{
    // VALUE ERROR
    if (mpz_sgn(op2) < 0)
        exit(EXIT_FAILURE);

    // Base case: if `op2` is 0, set `rop` to 1
    if (mpz_iszero(op2))
    {
        poly_set_ui(rop, 1);
        return;
    }

    // Base case: if `op2` is 1, set `rop` to `op1`
    if (mpz_isone(op2))
    {
        poly_set(rop, op1);
        return;
    }

    // Degenerate case: if `op2 >= 1` and `op1` is 0, set `rop` to 0
    if (poly_iszero(op1))
    {
        poly_set_ui(rop, 0);
        return;
    }

    // Degenerate case: if `op2 >= 1` and `op1` is 1, set `rop` to 1
    if (poly_isone(op1))
    {
        poly_set_ui(rop, 1);
        return;
    }

    // Otherwise, compute `op1 ** op2` with exponentiation by squaring
    POLY_GUARD_INIT(ptr, scratch, rop, rop == op1);

    /*************************************************************************
     * Computing P(x) ** p using exponentiation by squaring.                 *
     *                                                                       *
     * In basic terms, computing P(x) ** p is equivalent to computing either *
     *                                                                       *
     *   (P(x) ** 2) ** (p << 1),       if p is even, or                     *
     *   P(x) * (P(x) ** (p << 1)),     if p is odd.                         *
     *                                                                       *
     * Thus, we can reduce the problem of computing P(x) to the power of p,  *
     * with d = deg P(x), to the problem of computing (P(x) * P(x)) ** d/2   *
     * or P(x) * (P(x) ** (d - 1)/2). Either way, we compute a product of    *
     * polynomials of degree at most d, and a power of P(x) to floor(d/2).   *
     * Thus, it takes on the order of log(p) multiplications of polynomials  *
     * of degree at most d * p to compute P(x) ** p.                         *
     *************************************************************************/

    // Set `rop` to `op1**(2**b0) * op1**(2**b1) * ... * op1**(2**bn)`, where
    // `op2 == 2**b0 + 2**b1 + ... + 2**bn` is the binary expansion of `op2`.
    {
        poly_t base;
        poly_init_set(base, op1);

        mpz_t exp;
        mpz_init_set(exp, op2);
        
        poly_set_ui(ptr, 1);
        while (mpz_sgn(exp) > 0)
        {
            if (mpz_odd_p(exp))
                poly_mul(ptr, ptr, base);
            mpz_div_2exp(exp, exp, 1);
            if (mpz_cmp_ui(exp, 0) > 0)
                poly_mul(base, base, base);
        }

        mpz_clear(exp);
        poly_clear(base);
    }

    POLY_GUARD_CLEAR(rop, scratch);
}

// Sets `rop` to the power of `op1` by a non-negative integer `op2`. Faster
// than `poly_pow_z`.
void poly_pow_ui(poly_t rop, poly_t op1, unsigned long op2)
{
    // Base case: if `op2` is 0, set `rop` to 1
    if (ui_iszero(op2))
    {
        poly_set_ui(rop, 1);
        return;
    }

    // Base case: if `op2` is 1, set `rop` to `op1`
    if (ui_isone(op2))
    {
        poly_set(rop, op1);
        return;
    }

    // Degenerate case: if `op2 >= 1` and `op1` is 0, set `rop` to 0
    if (poly_iszero(op1))
    {
        poly_set_ui(rop, 0);
        return;
    }

    // Degenerate case: if `op2 >= 1` and `op1` is 1, set `rop` to 1
    if (poly_isone(op1))
    {
        poly_set_ui(rop, 1);
        return;
    }

    // Otherwise, compute `op1 ** op2` with exponentiation by squaring
    POLY_GUARD_INIT(ptr, scratch, rop, rop == op1);
    
    // Set `rop` to `op1**(2**b0) * op1**(2**b1) * ... * op1**(2**bn)`, where
    // `op2 == 2**b0 + 2**b1 + ... + 2**bn` is the binary expansion of `op2`.
    {
        poly_t base;
        poly_init_set(base, op1);

        poly_set_ui(ptr, 1);
        while (op2 > 0)
        {
            if (op2 & 1)
                poly_mul(ptr, ptr, base);
            op2 >>= 1;
            if (op2 > 0)
                poly_mul(base, base, base);
        }
        
        poly_clear(base);
    }

    POLY_GUARD_CLEAR(rop, scratch);
}

/* ************************************************************************* */