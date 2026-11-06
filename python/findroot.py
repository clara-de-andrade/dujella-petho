#!/usr/bin/python3

"""
This module contains a function findroot implementing the Newton-Raphson
algorithm for approximating the root of a function (in this case, a polynomial
function, where its coefficients are given) close to an initial value.
"""

from mpmath import mpf, fabs

def findroot(coeffs_in, init):
    """
    Computes an approximation x of the root of a polynomial function P, where

        P(x) = coeffs[0] + coeffs[1] * x + ... + coeffs[d] * x**d,

    with deg P = d, and x is close to init. In particular, P and its derivative
    dP are calculated using Horner's method.
    """

    # Initialise x
    x = mpf(init)
    # Horner's method requires the coefficients to be iterated in reverse order
    coeffs = list(reversed(coeffs_in))
    
    while True:
        # Computes P and its derivative dP using Horner's method. Note that:
        #
        #    P(x) = coeffs[0] + x * (coeffs[1] + x * (... + x * coeffs[d] ...))
        #
        # Therefore, P(x) is equal to P_d(x), where
        #
        #    P_0(x) = coeffs[d],
        #    P_(i + 1)(x) = coeffs[d - i] + x * P_i(x), 0 <= i < d.
        #
        # Moreover, by differentiating P_n(x), we have dP(x) = dP_d(x), where
        #
        #    dP_0(x) = 0,
        #    dP_(i + 1)(x) = P_i(x) + x * dP_i(x), 0 <= i < d.

        P  = coeffs[0]
        dP = 0
        for a in coeffs[1:]:
            dP = dP * x + P
            P  =  P * x + a

        step = P / dP
        tol = mp.eps * fabs(x)

        # Break loop if step is less than the tolerance for precision errors
        if fabs(step) <= tol:
            break

        # New value of x, given by the Newton-Raphson algorithm
        x_new = x - step
        
        # Break loop if x will not update to a new and distinct value
        if x_new == x:
            break

        # Update x to x_new
        x = x_new

    return x

