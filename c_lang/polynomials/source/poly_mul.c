#include "poly_ring.h"
#include "internal.h"

/* ************************************************************************* */

// Static helper function. Computes the product of `op1` and `op2` using the
// Cauchy product formula and sets it to `rop`.
static void poly_cauchy(poly_t rop, poly_t op1, poly_t op2)
{
    POLY_GUARD_INIT(ptr, scratch, rop, rop == op1 || rop == op2);

    // deg(op1 * op2) = deg(op1) + deg(op2)
    size_t deg = op1->deg + op2->deg;
    poly_set_deg(ptr, deg);

    /*************************************************************************
     * Computing P(x) * Q(x) by the Cauchy product formula.                  *
     *                                                                       *
     * If P and Q are given by                                               *
     *                                                                       *
     *   P(x) = a[0] + ... + x**m * a[m],                                    *
     *   Q(x) = b[0] + ... + x**n * b[n],                                    *
     *                                                                       *
     * then their product is given by                                        *
     *                                                                       *
     *   P(x) * Q)(x) = c[0] + ... + x**(m + n) * c[m + n],                  *
     *                                                                       *
     * where the coefficients c[0], ..., c[m + n] obey the formula           *
     *                                                                       *
     *   c[k] = a[0] * b[k] + a[1] * b[k - 1] + ... + a[k] * b[0]            *
     *                                                                       *
     * for every `0 <= k <= m + n`, setting `a[i] = 0`` and `b[j] = 0` when  *
     * `i > m` and `j > n`. When computing these coefficients in particular, *
     * we can prune the values of `a[i]` and `b[i]`, for `i` in the range    *
     * from `0` to `k` if we know them to always be 0 by this rule, i.e.     *
     * when k - i > n and i > m.                                             *
     *************************************************************************/

    // Set `prod(x)` to `op1(x) * op2(x)` by the Cauchy product formula.
    {
        // Temporary variable for coefficients of `prod`
        mpfr_t coeff;
        mpfr_init(coeff);
    
        // Set `prod->coeffs[k]` to `coeff` for each `0 <= k <= deg`
        for (size_t k = 0; k <= deg; k++)
        {
            // Initial value of `coeff` to 0
            mpfr_set_ui(coeff, 0, MPFR_RNDN);
    
            // Prune values in range `0 <= i <= k` where loop adds 0 to `coeff`
            size_t i_min = (k > op2->deg) ? (k - op2->deg) : 0;
            size_t i_max = min(k, op1->deg);

            for (size_t i = i_min; i <= i_max; i++)
                // Adds `op1->coeffs[i] * op2->coeffs[k - i]` to `coeff`
                mpfr_fma(coeff,
                    op1->coeffs[i], op2->coeffs[k - i], coeff, MPFR_RNDN);
            
            // Set `prod->coeffs[k]` to `coeff`
            mpfr_set(ptr->coeffs[k], coeff, MPFR_RNDN);
        }
    
        // Clear temporary variable
        mpfr_clear(coeff);
        poly_normalise(ptr);
    }

    POLY_GUARD_CLEAR(rop, scratch);
}

// Sets `rop` to the product of `op1` and `op2`. Computed for small degrees
// of `op1` and `op2` by the Cauchy product formula for formal power series,
// and by Karatsuba's algorithm for larger degrees.
void poly_mul(poly_t rop, poly_t op1, poly_t op2)
{
    // Degenerate case: if `op1` or `op2` are 0, set `rop` to 0
    if (poly_iszero(op1) || poly_iszero(op2))
    {
        poly_set_ui(rop, 0);
        return;
    }

    // Base case: if either `op1` or `op2` is constant, compute the product
    // of `op1` and `op2` by scalar multiplication (`poly_mul_fr`)
    if (op1->deg == 0)
    {
        // Degenerate case: if `op1` is 1, set `rop` to `op2`
        if (poly_isone(op1))
            poly_set(rop, op2);
        else
            poly_mul_fr(rop, op2, op1->coeffs[0]);
        return;
    }
    if (op2->deg == 0)
    {
        // Degenerate case: if `op2` is 1, set `rop` to `op1`
        if (poly_isone(op2))
            poly_set(rop, op1);
        else
            poly_mul_fr(rop, op1, op2->coeffs[0]);
        return;
    }

    // OVERFLOW ERROR
    if (op1->deg > SIZE_MAX - op2->deg)
        exit(EXIT_FAILURE);

    // Set `d` to the maximum of the lengths of `op1` and `op2`
    size_t d = max(op1->deg, op2->deg) + 1;

    // OVERFLOW ERROR
    if (d > (SIZE_MAX >> 1) + 1)
        exit(EXIT_FAILURE);

    // Otherwise, if `max(op1->deg, op2->deg)` is small, compute the product
    // of `op1` and `op2` using the Cauchy product formula (`poly_cauchy`)
    if (d <= POLY_MUL_CUTOFF)
    {
        poly_cauchy(rop, op1, op2);
        return;
    }

    // Otherwise, set `d` to the least power of 2 that is >= `d`
    d = 1ULL << (sizeof(unsigned long long) * 8 - __builtin_clzll(d - 1));

    /*************************************************************************
     * Computing P(x) * Q(x) by Karatsuba's algorithm.                       *
     *                                                                       *
     * Given P(x) and Q(x) polynomials of degree d, where d is a power of 2, *
     * let P0(x), P1(x) and Q0(x), Q1(x) be polynomials of degree d/2 where  * 
     *                                                                       *
     *   P(x) = P0(x) + x**(d/2) * P1(x),                                    *
     *   Q(x) = Q0(x) + x**(d/2) * Q1(x).                                    *
     *                                                                       *
     * Therefore, we have two identities involving P0, P1, Q0 and Q1:        *
     *                                                                       *
     *     P(x) * Q(x)                                                       *
     *   = (P0(x) + x**(d/2) * P1(x)) * (P0(x) + x**(d/2) * P1(x))           *
     *   =            (P0(x) * Q0(x))                                        *
     *   + x**(d/2) * ((P0(x) * Q1(x)) + (P1(x) * Q0(x)))                    *
     *   +     x**d * (P1(x) * Q1(x)),                                       *
     *                                                                       *
     *     (P0(x) + P1(x)) * (P0(x) + P1(x))                                 *
     *   = P0(x) * Q0(x)                                                     *
     *   + (P0(x) * Q1(x)) + (P1(x) * Q0(x))                                 *
     *   + P1(x) * Q1(x).                                                    *
     *                                                                       *
     * By the second, the term middle = (P0(x) * Q1(x)) + (P1(x) * Q0(x))    *
     * can be obtained by computing z1 - z0 - z2, where                      *
     *                                                                       *
     *   z0 = P0(x) * Q0(x),                                                 *
     *   z1 = (P0(x) + P1(x)) * (P0(x) + P1(x)),                             *
     *   z2 = P1(x) * Q1(x),                                                 *
     *                                                                       *
     * while P(x) * Q(x) is equal to z0 + x**(d/2) * middle + x**d * z2, so  *
     * we reduced the product of two polynomials of degree d to 3 products   *
     * and 2 sums of polynomials of degree d/2, as well as three polynomial  *
     * additions and two polynomial "shifts" (multiplying by x**p for some   *
     * non-negative integer p). Hence, it takes on the order of d**log2(3)   *
     * operations to compute P * Q, where log2(3) is approximately 1.585.    *
     *************************************************************************/

    size_t m = d >> 1;

    // Declare P0, P1, Q0, Q1 such as above 
    poly_t P0, P1, Q0, Q1;
    poly_inits(P0, P1, Q0, Q1, (struct poly_t*) NULL);
    {
        poly_realloc(P0, m);
        poly_realloc(P1, m);
        for (size_t i = 0; i <= op1->deg; i++)
            if (i < m)
                poly_set_coeff_fr(P0, op1->coeffs[i], i);
            else
                poly_set_coeff_fr(P1, op1->coeffs[i], i - m);
    
        poly_realloc(Q0, m);
        poly_realloc(Q1, m);
        for (size_t i = 0; i <= op2->deg; i++)
            if (i < m)
                poly_set_coeff_fr(Q0, op2->coeffs[i], i);
            else
                poly_set_coeff_fr(Q1, op2->coeffs[i], i - m);
    }

    // Declare z0, z1, z2 and middle such as above
    poly_t z0, z1, z2, middle;
    poly_inits(z0, z1, z2, middle, (struct poly_t*) NULL);
    {
        // z0 = P0 * Q0
        poly_mul(z0, P0, Q0);
        // z2 = P0 * Q0
        poly_mul(z2, P1, Q1);

        // z1 == (P0 * P1) + (Q0 * Q1)
        poly_add(z1, P0, P1);
        poly_add(middle, Q0, Q1);  // Using `middle` as temporary variable
        poly_mul(z1, z1, middle);  // Using `middle` as temporary variable
    
        // middle = z1 - z0 - z2
        poly_add(middle, z0, z2);
        poly_sub(middle, z1, middle);
    }
    poly_clears(P0, P1, Q0, Q1, (struct poly_t*) NULL);

    // Set `rop` to `z0 + (middle << m) + (z2 << d)`.
    {
        size_t deg = op1->deg + op2->deg;
        poly_set_deg(rop, deg);

        for (size_t i = 0; i <= deg; i++)
            if (i <= z0->deg)
                mpfr_set(rop->coeffs[i], z0->coeffs[i], MPFR_RNDN);
            else
                mpfr_set_ui(rop->coeffs[i], 0, MPFR_RNDN);

        for (size_t i = 0; i <= middle->deg; i++)
            mpfr_add(rop->coeffs[i + m],
                rop->coeffs[i + m], middle->coeffs[i], MPFR_RNDN);

        for (size_t i = 0; i <= z2->deg; i++)
            mpfr_add(rop->coeffs[i + d],
                rop->coeffs[i + d], z2->coeffs[i], MPFR_RNDN);
        
        poly_normalise(rop);
    }
    poly_clears(z0, z1, z2, middle, (struct poly_t*) NULL);
}

/* ************************************************************************* */

// Sets `rop` to the scalar product of `op1` by `op2`.
void poly_mul_fr(poly_t rop, poly_t op1, mpfr_t op2)
{
    if (poly_iszero(op1) || mpfr_zero_p(op2))
    {
        poly_set_ui(rop, 0);
        return;
    }

    POLY_SCALAR_GUARD_INIT(ptr, scratch, rop, op2);

    poly_set_deg(rop, op1->deg);
    for (size_t i = 0; i <= op1->deg; i++)
        mpfr_mul(rop->coeffs[i], op1->coeffs[i], ptr, MPFR_RNDN);
    poly_normalise(rop);

    POLY_SCALAR_GUARD_CLEAR(scratch);
}

// Sets `rop` to the scalar product of `op1` by `op2`.
#define DEFINE_POLY_MUL(SUFFIX, TYPE, MUL_METHOD, IS_ZERO)                  \
    void poly_mul_##SUFFIX(poly_t rop, poly_t op1, TYPE op2)                \
    {                                                                       \
        if (poly_iszero(op1) || IS_ZERO(op2))                        \
        {                                                                   \
            poly_set_ui(rop, 0);                                            \
            return;                                                         \
        }                                                                   \
                                                                            \
        poly_set_deg(rop, op1->deg);                                        \
        for (size_t i = 0; i <= op1->deg; i++)                              \
            MUL_METHOD(rop->coeffs[i], op1->coeffs[i], op2, MPFR_RNDN);     \
        poly_normalise(rop);                                                \
    }

DEFINE_POLY_MUL(z,      mpz_t,          mpfr_mul_z,     mpz_iszero)
DEFINE_POLY_MUL(ui,     unsigned long,  mpfr_mul_ui,    ui_iszero)
DEFINE_POLY_MUL(si,     long,           mpfr_mul_si,    si_iszero)

/* ************************************************************************* */

// Shifts `op1` by `op2` places and sets to `rop`. Equivalent to setting
// `rop(x)` to the product of `op1(x)` by `x ** op2`.
void poly_shift(poly_t rop, poly_t op1, size_t op2)
{
    // Degenerate case: if `op1` is 0, set `rop` to 0
    if (poly_iszero(op1))
    {
        poly_set_ui(rop, 0);
        return;
    }

    size_t old_deg = op1->deg;
    // OVERFLOW ERROR
    if (op2 > SIZE_MAX - old_deg)
        exit(EXIT_FAILURE);
    // deg(op1 << op2) = deg(op1) + op2
    poly_set_deg(rop, old_deg + op2);

    // Shift coefficients of `op1` by `op2` and set to `rop`
    for (size_t i = old_deg; i != (size_t) -1; i--)
        mpfr_set(rop->coeffs[i + op2], op1->coeffs[i], MPFR_RNDN);

    // Zero out the coefficients of `rop` less than `op2`
    for (size_t i = 0; i < op2; i++)
        mpfr_set_ui(rop->coeffs[i], 0, MPFR_RNDN);
}

/* ************************************************************************* */