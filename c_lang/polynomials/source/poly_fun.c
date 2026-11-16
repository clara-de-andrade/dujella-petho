#include "poly_ring.h"
#include "internal.h"

/* ************************************************************************* */

// Static helper function. Recursor for `poly_eval`.
static void poly_eval_rec(poly_t rop, poly_t op1, poly_t* op2_pow, size_t k)
{
    POLY_GUARD_INIT(ptr, scratch, rop, rop == op1);

    // Degenerate case: if `op1` is constant, return a deep copy of `op1`
    if (op1->deg == 0)
    {
        poly_set_fr(ptr, op1->coeffs[0]);
        poly_normalise(ptr);
        goto end;
    }

    // Otherwise, run recursive algorithm for computing `op1 ** op2`
    size_t m = (size_t) 1 << (k - 1);
    
    /*************************************************************************
     * Divide and conquer algorithm for computing P(Q(x)).                   *
     *                                                                       *
     * To compute the composition of the polynomials P(x) and Q(x) of degree *
     * d and n, where both are powers of 2, is equivalent to evaluating P(t) *
     * at t = Q(x). First, we must find two polynomials P0(t) and P1(t) of   *
     * degree d/2 such that                                                  *
     *                                                                       *
     *   P(t) = P0(t) + (t**(d/2)) * P1(t).                                  *
     *                                                                       *
     * Thus, we have reduced the problem of evaluating a degree d polynomial *
     * to that of evaluating two degree d/2 polynomials. In addition, we've  *
     * raised t = Q(x) to the power of d/2, which involves on the order of   *
     * log(d/2) polynomial multiplications of degree <= (d/2) * n, and then  *
     * we multiply the result by P1(t). By Karatsuba's algorithm to compute  *
     *  products of polynomials, it takes on the order of (n * d)**log2(3)   *
     * operations to compute P(t) at t = Q(x), and log2(3) is approximately  *
     * equal to 1.585.                                                       *
     *************************************************************************/
    
    // Declare `P0`, `P1` such that `op1(x) = P0(x) + (x ** (d/2)) * P1(x)`
    poly_t P0, P1;
    poly_inits(P0, P1, (struct poly_t*) NULL);
    {
        poly_realloc(P0, m);
        poly_realloc(P1, m);
        for (size_t i = 0; i <= op1->deg; i++)
            if (i < m)
                poly_set_coeff_fr(P0, op1->coeffs[i], i);
            else
                poly_set_coeff_fr(P1, op1->coeffs[i], i - m);
    }

    // Evaluate P0(x) + (x ** (d/2)) * P1(x) at `x = op2` and set to `rop`.
    {
        poly_t P0_eval, P1_eval;
        poly_inits(P0_eval, P1_eval, (struct poly_t*) NULL);

        poly_eval_rec(P0_eval, P0, op2_pow, k - 1);
        poly_eval_rec(P1_eval, P1, op2_pow, k - 1);
    
        poly_mul(P1_eval, P1_eval, op2_pow[k - 1]);
        poly_add(ptr, P0_eval, P1_eval);

        poly_clears(P0_eval, P1_eval, (struct poly_t*) NULL);
    }
    poly_clears(P0, P1, (struct poly_t*) NULL);

    end:
        POLY_GUARD_CLEAR(rop, scratch);
}

// Sets `rop` to the composition of `op1` and `op2`.
void poly_eval(poly_t rop, poly_t op1, poly_t op2)
{
    POLY_GUARD_INIT(ptr, scratch, rop, rop == op1 || rop == op2);
    
    // Degenerate case: if `op1` is constant, return a deep copy of `op1`
    if (op1->deg == 0)
    {
        poly_set_fr(ptr, op1->coeffs[0]);
        poly_normalise(ptr);

        goto end;
    }

    // Base case: if `op2` is a constant, evaluate `op1` at `op2->coeff[0]`
    if (op2->deg == 0)
    {
        mpfr_t val;
        mpfr_init(val);

        poly_eval_fr(val, op1, op2->coeffs[0]);
        poly_set_fr(ptr, val);

        mpfr_clear(val);
        goto end;
    }

    // OVERFLOW ERROR
    if (op1->deg > SIZE_MAX / op2->deg)
        exit(EXIT_FAILURE);

    // Otherwise, set `d` to the least power of 2 that is >= `op1->deg + 1`
    size_t d = op1->deg + 1;
    // OVERFLOW ERROR
    if (d > (SIZE_MAX >> 1) + 1)
        exit(EXIT_FAILURE);
    // Declare `k` such that `d = 2**k`.
    size_t k = sizeof(unsigned long long) * 8 - __builtin_clzll(d - 1);

    // Temporary array of variables for powers of op_2
    poly_t* op2_pow = malloc(k * sizeof(poly_t));
    // MEMORY ERROR
    if (op2_pow == NULL)
        exit(EXIT_FAILURE);

    // Set `op2_pow[i]` to `op2**(2**i)` for each `0 <= i < k`
    poly_init_set(op2_pow[0], op2);
    for (size_t j = 1; j < k; j++)
    {
        poly_init(op2_pow[j]);
        poly_mul(op2_pow[j], op2_pow[j - 1], op2_pow[j - 1]);
    }

    // Run recursor with precomputed powers of `op2`
    poly_eval_rec(ptr, op1, op2_pow, k);

    // Clear temporary variables and free array
    for (size_t j = 0; j < k; j++)
        poly_clear(op2_pow[j]);
    free(op2_pow);

    end:
        POLY_GUARD_CLEAR(rop, scratch);
}

/* ************************************************************************* */

// Evaluates `op1(x)` at `x = op2` using Horner' method and sets to `rop.
void poly_eval_fr(mpfr_t rop, poly_t op1, mpfr_t op2)
{
    /*************************************************************************
     * Horner's method for evaluating a polynomial P(x) at x = t.            *
     *                                                                       *
     * If P(x) = a(0) + ... + x**n * a(n) is a polynomial, then we have that *
     *                                                                       *
     *   P(x) = a(0) + x * (a(1) + x * (... (a(n - 1) + x * a(m)) ...));     *
     *                                                                       *
     * that is P(x) = P(n)(x), where P(0)(x), ..., P(n)(x) are given by      *
     *                                                                       *
     *   P(0)(x) = a(n),                                                     *
     *   P(i)(x) = x * P(i - 1)(x) + a(n - i), 0 <= i <= n,                  *
     *                                                                       *
     * Thus, instead of computing n - 1 powers t**2, ..., t**n of t, as well *
     * as n - 1 products t * a(1), ..., t**n * a(n) and summing the products *
     * to obtain P(x) at x = t, we instead compute n products by t, namely   *
     * t * P(0)(t), ..., t * P(n - 1)(t) and n sums of each product by the   *
     * coefficients a(1), ..., a(n), respectively. In total, it takes on the *
     * order of n operations (sums and scalar products) to evaluate P(t).    *
     *************************************************************************/

    mpfr_t acc;
    mpfr_init_set_ui(acc, 0, MPFR_RNDN);

    for (size_t i = 0; i <= op1->deg; i++)
        mpfr_fma(acc, op2, acc, op1->coeffs[op1->deg - i], MPFR_RNDN);
    mpfr_set(rop, acc, MPFR_RNDN);
    mpfr_clear(acc);
}

// Evaluates `op1(x)` at `x = op2` using Horner' method and sets to `rop.
#define DEFINE_POLY_EVAL(SUFFIX, TYPE, INIT_SET_METHOD)                         \
    void poly_eval_##SUFFIX(mpfr_t rop, poly_t op1, TYPE op2)                   \
    {                                                                           \
        mpfr_t acc;                                                             \
        mpfr_init_set_ui(acc, 0, MPFR_RNDN);                                    \
                                                                                \
        mpfr_t op2_fr;                                                          \
        INIT_SET_METHOD(op2_fr, op2, MPFR_RNDN);                                \
                                                                                \
        for (size_t i = 0; i <= op1->deg; i++)                                  \
            mpfr_fma(acc, op2_fr, acc, op1->coeffs[op1->deg - i], MPFR_RNDN);   \
        mpfr_set(rop, acc, MPFR_RNDN);                                          \
                                                                                \
        mpfr_clears(acc, op2_fr, (mpfr_ptr) NULL);                              \
    }


DEFINE_POLY_EVAL(z,     mpz_t,          mpfr_init_set_z)
DEFINE_POLY_EVAL(ui,    unsigned long,  mpfr_init_set_ui)
DEFINE_POLY_EVAL(si,    long,           mpfr_init_set_si)

/* ************************************************************************* */

// Sets `rop` to derivative of `op`, computed by symbolic differentiation.
void poly_d(poly_t rop, poly_t op)
{
    // Degenerate case: if `op` is 0, set `rop` to 0
    if (op->deg == 0)
    {
        poly_set_ui(rop, 0);
        return;
    }

    size_t orig_deg = op->deg;
    // Allocates at least `op->deg` coefficients to `rop`
    poly_realloc(rop, orig_deg);

    // Computes coefficients of derivative of `op` and sets to `rop`
    for (size_t i = 1; i <= orig_deg; i++)
        mpfr_mul_ui(rop->coeffs[i - 1], op->coeffs[i], i, MPFR_RNDN);

    poly_set_deg(rop, orig_deg - 1);
    poly_normalise(rop);
}

/* ************************************************************************* */