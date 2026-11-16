#include "internal.h"

// TODO: poly_mullow

/* ************************************************************************* */

// Static helper function. Computes the product of `op1` and `op2` modulo
// `x**d` using the Cauchy product formula and sets it to `rop`.
static void poly_cauchy(poly_t rop, poly_t op1, poly_t op2, size_t d)
{
    POLY_GUARD_INIT(ptr, scratch, rop, rop == op1 || rop == op2);

    // TODO: document
    size_t deg1 = min(op1->deg, d - 1);
    size_t deg2 = min(op2->deg, d - 1);

    size_t deg = min(deg1 + deg2, d - 1);
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
        // Temporary variable
        mpfr_t acc;
        mpfr_init(acc);
    
        // Set `prod->coeffs[k]` to `acc` for each `0 <= k <= deg`
        for (size_t k = 0; k <= deg; k++)
        {
            // Initial value of `acc` to 0
            mpfr_set_ui(acc, 0, MPFR_RNDN);
    
            // Prune values in range `0 <= i <= k` where loop adds 0 to `acc`
            size_t i_min = (k > deg2) ? (k - deg2) : 0;
            size_t i_max = min(k, deg1);

            for (size_t i = i_min; i <= i_max; i++)
                // Adds `op1->coeffs[i] * op2->coeffs[k - i]` to `acc`
                mpfr_fma(acc,
                    op1->coeffs[i], op2->coeffs[k - i], acc, MPFR_RNDN);
            
            // Set `prod->coeffs[k]` to `acc`
            mpfr_set(ptr->coeffs[k], acc, MPFR_RNDN);
        }
    
        // Clear temporary variable
        mpfr_clear(acc);
        poly_normalise(ptr);
    }

    POLY_GUARD_CLEAR(rop, scratch);
}

// Sets `rop` to the product of `op1` and `op2` modulo `x**d`.
void poly_mullow(poly_t rop, poly_t op1, poly_t op2, size_t d)
{
    // Degenerate case: if `op1` or `op2` are 0, set `rop` to 0
    if (d == 0 || poly_zero(op1) || poly_zero(op2))
    {
        poly_set_ui(rop, 0);
        return;
    }

    // TODO: document
    if (d == 1)
    {
        mpfr_t coeff;
        mpfr_init(coeff);

        mpfr_mul(coeff, op1->coeffs[0], op2->coeffs[0], MPFR_RNDN);
        poly_set_fr(rop, coeff);

        mpfr_clear(coeff);
        return;
    }

    // Base case: if either `op1` or `op2` is constant, compute the product
    // of `op1` and `op2` by scalar multiplication (`poly_mul_fr`)
    if (op1->deg == 0)
    {
        // Degenerate case: if `op1` is 1, set `rop` to `op2`
        if (poly_unit(op1))
            poly_set(rop, op2);
        else
            poly_mul_fr(rop, op2, op1->coeffs[0]);

        if (rop->deg >= d)
        {
            poly_set_deg(rop, d - 1);
            poly_normalise(rop);
        }
        return;
    }
    if (op2->deg == 0)
    {
        // Degenerate case: if `op2` is 1, set `rop` to `op1`
        if (poly_unit(op2))
            poly_set(rop, op1);
        else
            poly_mul_fr(rop, op1, op2->coeffs[0]);

        if (rop->deg >= d)
        {
            poly_set_deg(rop, d - 1);
            poly_normalise(rop);
        }
        return;
    }

    // Otherwise, if `d` is small enough, compute the product of `op1` and 
    // op2` using the Cauchy product formula (`poly_cauchy`)
    if (d <= POLY_MUL_CUTOFF)
    {
        poly_cauchy(rop, op1, op2, d);
        return;
    }

    size_t deg1 = min(op1->deg, d - 1);
    size_t deg2 = min(op2->deg, d - 1);

    // Evaluates to true if computing `op1 * op2 mod x**d`, and to false if
    // computing `op1 * op2` without truncation. 
    int full_product =
        (op1->deg <= SIZE_MAX - op2->deg && d - 1 >= op1->deg + op2->deg);

    size_t m;
    if (full_product)
    {   // Otherwise, if `d` is large enough, set `d` to the least power of 2
        // that is >= the greatest of the lengths of `op1` and `op2`.
        d = max(op1->deg, op2->deg) + 1;
        // OVERFLOW ERROR
        if (d > (SIZE_MAX >> 1) + 1)
            exit(EXIT_FAILURE);        
        d = 1ULL << (sizeof(unsigned long long) * 8 - __builtin_clzll(d - 1));

        // ...and set `m` to `d/2`.
        m = d >> 1;
    }
    else
        // Otherwise, if `d` is large enough, set `m` to `ceil(d/2)`.
        m = (d >> 1) + (d & 1);

    // TODO: update documentation

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

    // Declare P0, P1, Q0, Q1 such as above 
    poly_t P0, P1, Q0, Q1;
    poly_inits(P0, P1, Q0, Q1, (struct poly_t*) NULL);
    {
        for (size_t i = 0; i <= deg1; i++)
            if (i < m)
                poly_set_coeff_fr(P0, op1->coeffs[i], i);
            else
                poly_set_coeff_fr(P1, op1->coeffs[i], i - m);

        for (size_t i = 0; i <= deg2; i++)
            if (i < m)
                poly_set_coeff_fr(Q0, op2->coeffs[i], i);
            else
                poly_set_coeff_fr(Q1, op2->coeffs[i], i - m);
    }

    if (full_product)
    {
        // Declare z0, z1, z2 and middle such as above
        poly_t z0, z1, z2, middle;
        poly_inits(z0, z1, z2, middle, (struct poly_t*) NULL);

        // z0 = P0 * Q0
        poly_mul(z0, P0, Q0);

        // z1 == (P0 * P1) + (Q0 * Q1)
        poly_add(z1, P0, P1);
        poly_add(middle, Q0, Q1);  // Using `middle` as temporary variable
        poly_mul(z1, z1, middle);  // Using `middle` as temporary variable
        
        // z2 = P1 * Q1
        poly_mul(z2, P1, Q1);
        
        // middle = z1 - z0 - z2
        poly_add(middle, z0, z2);
        poly_sub(middle, z1, middle);

        size_t deg = deg1 + deg2;
        poly_set_deg(rop, deg);

        // Set `rop` to `z0 + (middle << m) + (z2 << d)`.
        for (size_t i = 0; i <= deg; i++)
            if (i <= z0->deg)
                mpfr_set(rop->coeffs[i], z0->coeffs[i], MPFR_RNDN);
            else
                mpfr_set_ui(rop->coeffs[i], 0, MPFR_RNDN);

        if (!poly_zero(middle))
            for (size_t i = 0; i <= middle->deg; i++)
                mpfr_add(rop->coeffs[i + m],
                    rop->coeffs[i + m], middle->coeffs[i], MPFR_RNDN);

        if (!poly_zero(z2))
            for (size_t i = 0; i <= z2->deg; i++)
                mpfr_add(rop->coeffs[i + d],
                    rop->coeffs[i + d], z2->coeffs[i], MPFR_RNDN);

        poly_clears(z0, z1, z2, middle, (struct poly_t*) NULL);
    }
    else
    {
        // TODO: document
        poly_t z0, middle, temp;
        poly_inits(z0, middle, temp, (struct poly_t*) NULL);
        
        // z0 = P0 * Q0
        poly_mul(z0, P0, Q0);

        // middle = (P1 * Q0) + (P0 * Q1)
        poly_mullow(temp, P0, Q1, d - m);
        poly_mullow(middle, P1, Q0, d - m);
        poly_add(middle, temp, middle);

        size_t deg = deg1 + deg2;
        poly_set_deg(rop, deg);

        for (size_t i = 0; i <= deg; i++)
            if (i <= z0->deg)
                mpfr_set(rop->coeffs[i], z0->coeffs[i], MPFR_RNDN);
            else
                mpfr_set_ui(rop->coeffs[i], 0, MPFR_RNDN);

        if (!poly_zero(middle))
        {
            size_t i_max = min(middle->deg, deg - m);
            for (size_t i = 0; i <= i_max; i++)
                mpfr_add(rop->coeffs[i + m],
                    rop->coeffs[i + m], middle->coeffs[i], MPFR_RNDN);
        }

        poly_clears(z0, middle, temp, (struct poly_t*) NULL);
    }

    poly_normalise(rop);

    poly_clears(P0, P1, Q0, Q1, (struct poly_t*) NULL);
}

/* ************************************************************************* */

// Sets `rop` to the scalar product of `op1` by `op2`.
void poly_mul_fr(poly_t rop, poly_t op1, mpfr_t op2)
{
    if (poly_zero(op1) || mpfr_zero_p(op2))
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
        if (poly_zero(op1) || IS_ZERO(op2))                        \
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

DEFINE_POLY_MUL(z,      mpz_t,          mpfr_mul_z,     mpz_zero)
DEFINE_POLY_MUL(ui,     unsigned long,  mpfr_mul_ui,    ui_zero)
DEFINE_POLY_MUL(si,     long,           mpfr_mul_si,    si_zero)

/* ************************************************************************* */

// Shifts `op1` by `op2` places and sets to `rop`. Equivalent to setting
// `rop(x)` to the product of `op1(x)` by `x ** op2`.
void poly_shift(poly_t rop, poly_t op1, size_t op2)
{
    // Degenerate case: if `op1` is 0, set `rop` to 0
    if (poly_zero(op1))
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