#include <stdio.h>
#include "polynom.h"

/* ************************************************************************* */

// Creates polynomial P(x) = coeffs[0] + ... + coeffs[deg - 1] * x**(deg - 1).
poly_t* poly_init(size_t deg, mpz_t* coeffs)
{
    poly_t* poly = malloc(sizeof(poly_t));

    poly->deg = deg;
    poly->coeffs = calloc(deg, sizeof(mpz_t));

    for (size_t i = 0; i < deg; i++)
        mpz_init_set(poly->coeffs[i], coeffs[i]);

    return poly;
}

void poly_clear(poly_t* P)
{
    // TODO
}

/* ************************************************************************* */

// Computes the sum of two polynomials A(t) = a_0 + ... + a_m * t**m and
// B = b_0 + ... + b_n * t**n, where m, n >= 0, by pointwise addition.
poly_t* poly_add(poly_t* A, poly_t* B)
{
    size_t deg = max(a->deg, b->deg);
    mpz_t coeffs[deg];

    for (size_t i = 0; i < deg; i++)
    {
        mpz_init_set(coeffs[i], 0);
        if (i < A->deg)
            mpz_add(coeffs[i], coeffs[i], a->coeffs[i]);
        if (i < B->deg)
            mpz_add(coeffs[i], coeffs[i], b->coeffs[i]);
    }

    poly_t* res = poly_init(deg, coeffs);
    return res;
}

#define max(a, b) ((a) >= (b) ? (a) : (b))

// Computes the product of two polynomials A(t) = a_0 + ... + a_m * t**m and
// B = b_0 + ... + b_n * t**n, where m, n >= 0, using Karatsuba's algorithm.
poly_t* poly_mul(poly_t* A, poly_t* B)
{
    size_t d = 1 << mpz_sizeinbase(max(A->deg, B->deg), 2);
    if (d == 1)
    {
        mpz_t coeffs[1];
        mpz_init(coeffs[0]);
        mpz_mul(coeffs[0], A->coeffs[0], B->coeffs[0]);

        poly_t* res = poly_init(1, coeffs);
        return res;
    }

    size_t m = d >> 1;
    mpz_t A0_coeffs[m];
    mpz_t A1_coeffs[m];
    mpz_t B0_coeffs[m];
    mpz_t B1_coeffs[m];

    for (size_t i = 0; i < m; i++)
    {
        mpz_init(A0_coeffs[i]);
        if (i < A->deg)
            mpz_set(A0_coeffs[i], A->coeffs[i]);
        else
            mpz_set_ui(A0_coeffs[i], 0);
        
        mpz_init(B0_coeffs[i]);
        if (i < B->deg)
            mpz_set(B0_coeffs[i], B->coeffs[i]);
        else
            mpz_set_ui(B0_coeffs[i], 0);
    }

    poly_t* A0 = poly_init(m, A0_coeffs);
    poly_t* B0 = poly_init(m, B0_coeffs);

    for (size_t i = m; i < d; i++)
    {
        mpz_init(A1_coeffs[i]);
        if (i < A->deg)
            mpz_set(A1_coeffs[i], A->coeffs[i - m]);
        else
            mpz_set_ui(A1_coeffs[i], 0);
        
        mpz_init(B1_coeffs[i]);
        if (i < B->deg)
            mpz_set(B1_coeffs[i], B->coeffs[i - m]);
        else
            mpz_set_ui(B1_coeffs[i], 0);
    }

    poly_t* A1 = poly_init(m, A1_coeffs);
    poly_t* B1 = poly_init(m, B1_coeffs);


    poly_t* z0 = polymul(A0, B0);
    poly_t* z1 = polymul(polyadd(A0, A1), polyadd(B0, B1));
    poly_t* z2 = polymul(A1, B1);

    mpz_t z3_coeffs[z1->deg];
    for (size_t i = 0; i < z1->deg; i++)
    {
        mpz_init_set(z3_coeffs[i], 0);
        mpz_sub(z3_coeffs[i], z1->coeffs[i], z0_coeffs[i]);
        mpz_sub(z3_coeffs[i], z3->coeffs[i], z2_coeffs[i]);
    }
    poly_t* z3 = poly_init(z1->deg, z3_coeffs);

    mpz_t coeffs[d << 1];
    for (size_t i = 0; i < z0->deg; i++)
        mpz_init_set(coeffs[i], z0->coeffs[i]);
    for (size_t i = 0; i < z3->deg; i++)
        mpz_init_set(coeffs[m + i], z3->coeffs[i]);
    for (size_t i = 0; i < z2->deg; i++)
        mpz_init_set(coeffs[d + i], z2->coeffs[i]);

    poly_t* res = poly_init(deg, coeffs);
    return res;
}

/* ************************************************************************* */

poly_t* poly_diff(poly_t* P)
{
    // TODO
}

/* ************************************************************************* */

// Evaluates polynomial P at real value x.
void poly_eval(mpfr_t op, poly_t* P, mpfr_t x)
{
    // TODO
}

/* ************************************************************************* */
