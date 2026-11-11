#include <stdio.h>
#include "polynom.h"

/* ************************************************************************* */

// C Type of polynomials with integer coefficients
struct _poly_t {
    size_t deg;
    mpz_t* coeffs;
};

// Creates a new polynomial coeffs[0] + ... + coeffs[deg] * x**(deg)
poly_t* poly_new(size_t deg, mpz_t* coeffs)
{
    poly_t* poly = malloc(sizeof(poly_t));

    if (poly == (poly_t*) NULL)
        goto end;

    poly->deg = deg;

    if (coeffs == (mpz_t*) NULL)
    {
        coeffs = calloc(deg + 1, sizeof(mpz_t));
        for (size_t i = 0; i <= deg; i++)
        {
            mpz_init(coeffs[i]);
            mpz_set_ui(coeffs[i], 0);
        }
    }

    poly->coeffs = calloc(deg, sizeof(mpz_t));

    if (poly->coeffs == (mpz_t*) NULL || coeffs == (mpz_t*) NULL)
        goto end;

    for (size_t i = 0; i < deg; i++)
        mpz_init_set(poly->coeffs[i], coeffs[i]);

    end:
        return poly;
}

// Frees self
void poly_free(poly_t* self)
{
    if (poly == (poly_t*) NULL)
        return;

    if (poly->coeffs == (mpz_t*) NULL)
        goto end;
    
    for (size_t i = 0; i < poly->deg; i++)
        mpz_clear(poly->coeffs[i]);

    end:
        free(poly);
}

/* ************************************************************************* */

// Getter for deg
size_t poly_get_deg(poly_t* self);
// Getter for coeffs
mpz_t* poly_get_coeffs(poly_t* self);
// Getter for coeffs[key]
void poly_get_coeff(poly_t* self, size_t key, mpz_t op);

// Setter for deg
void poly_set_deg(poly_t* self, size_t value);
// Setter for coeffs
void poly_set_coeffs(poly_t* self, mpz_t* value);
// Setter for coeffs[key]
void poly_set_coeff(poly_t* self, size_t key, mpz_t value);

/* ************************************************************************* */

// String representation of self
char* poly_str(poly_t* self);
// Number of coefficients in self
size_t poly_len(poly_t* self);

// Evaluate self(x) at x = other
int64_t poly_eval_si(poly_t* self, int64_t other);
uint64_t poly_eval_ui(poly_t* self, uint64_t other);
void poly_eval_z(poly_t* self, mpz_t other, mpz_t op);
void poly_eval_fr(poly_t* self, mpfr_t other, mpfr_t op);
poly_t* polt_eval(poly_t* self, poly_t* other);

/* ************************************************************************* */

// Computes self(x) + other(x) by pointwise addition
poly_t* poly_add(poly_t* self, poly_t* other)
{
    size_t deg = max(A->deg, B->deg);
    mpz_t coeffs[deg];

    for (size_t i = 0; i < deg; i++)
    {
        mpz_init_set(coeffs[i], 0);
        if (i < A->deg)
            mpz_add(coeffs[i], coeffs[i], a->coeffs[i]);
        if (i < B->deg)
            mpz_add(coeffs[i], coeffs[i], b->coeffs[i]);
    }

    poly_t* res = poly_new(deg, coeffs);
    return res;
}

// Computes self(x) - other(x) by pointwise addition
poly_t* poly_sub(poly_t* self, poly_t* other);

// Computes other * self(x) by scalar multiplication
poly_t* poly_mul_si(poly_t* self, int64_t other);
poly_t* poly_mul_ui(poly_t* self, uint64_t other);
poly_t* poly_mul_z(poly_t* self, mpz other);
poly_t* poly_mul_fr(poly_t* self, mpfr_t other);

#define max(a, b) ((a) >= (b) ? (a) : (b))

// Computes self(x) * other(x) using Karatsuba's algorithm
poly_t* poly_mul(poly_t* self, poly_t* other)
{
    size_t d = 1 << mpz_sizeinbase(max(A->deg, B->deg), 2);
    if (d == 1)
    {
        mpz_t coeffs[1];
        mpz_init(coeffs[0]);
        mpz_mul(coeffs[0], A->coeffs[0], B->coeffs[0]);

        poly_t* res = poly_new(1, coeffs);
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

    poly_t* A0 = poly_new(m, A0_coeffs);
    poly_t* B0 = poly_new(m, B0_coeffs);

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

    poly_t* A1 = poly_new(m, A1_coeffs);
    poly_t* B1 = poly_new(m, B1_coeffs);


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
    poly_t* z3 = poly_new(z1->deg, z3_coeffs);

    mpz_t coeffs[d << 1];
    for (size_t i = 0; i < z0->deg; i++)
        mpz_init_set(coeffs[i], z0->coeffs[i]);
    for (size_t i = 0; i < z3->deg; i++)
        mpz_init_set(coeffs[m + i], z3->coeffs[i]);
    for (size_t i = 0; i < z2->deg; i++)
        mpz_init_set(coeffs[d + i], z2->coeffs[i]);

    poly_t* res = poly_new(deg, coeffs);
    return res;
}

// Computes self(x) ** exp using exponentiation by squaring
poly_t* poly_pow_ui(poly_t* self, uint64_t exp);
poly_t* poly_pow_z(poly_t* self, mpz_t exp);

// Computes -self(x) = (-1) * self(x)
poly_t* poly_neg(poly_t* self);

/* ************************************************************************* */

// Differentiates of polynomial P(t) = a_0 + ... + a_n * t**n, where n >= 0.
poly_t* poly_diff(poly_t* P);

/* ************************************************************************* */
