#!/usr/bin/python3

"""
This module contains a function reqseq which computes the n-th term of a linear
recurrent sequence F_n, given by its coefficients and its initial values, where

    F_n = coeffs[0] * F_(n - 1) + ... + coeffs[k - 1] F_(n - k),       (n >= k)
    F_0 = inits[0], ..., F_(k - 1) = inits[k - 1]

The implementation of reqseq uses the generating function of F_n, given by

    f(q, k)(t) = P(t), Q(t),

for all t in its radius of convergence, where

    P(t) = sum(
        0 <= n <= k - 1, 
        (inits[n] - sum(1 <= i <= n, coeffs[i] * inits[n - i])) * x**n
    ),
    Q(t) = 1 - sum(0 <= n <= k - 1, coeffs[n] * x**n).
"""

def reqseq(k, coeffs, inits, n):
    """
    Computes the n-th term of a linear recurrent sequence, as the n-th
    coefficient of the formal power series for its generating function,
    using the Bostan-Mori algorithm. (Algorithm 1 of arXiv:2008.08822v1.)
    """

    if n >= 0 and n < k:
        return inits[n]
    
    P = sum(
        (inits[n] - sum(
            coeffs[i] * inits[n - i]
            for i in range(1, n + 1)))
        for n in range(k))
    Q = [1] + [-c for c in coeffs]

    while n > 0:
        # Q_minus(t) = Q(-t)
        Q_minus = [Q[i] if i % 2 == 0 else -Q[i] for i in range(k + 1)]
        # U(t) = P(t) * Q(-t)
        U = polymul(P, Q_minus)
        # V(t) = Q(t) * Q(-t)
        V = polymul(Q, Q_minus)

        # Note that, since P(t) / Q(t) = U(t) / V(t), if we let
        #
        #     U(t) = A(t**2) + t * B(t**2),
        #     m = floor(k / 2)
        #
        # then the n-th coefficient of P(t) / Q(t) is equal to the m-th
        # coefficient of A(t)/V(t) if n is even, and equal to the m-th
        # coefficient of B(t)/V(t) if n is odd.

        if n % 2 == 0:
            P = U[0::2][:k]
        else:
            P = U[1::2][:k]
        Q = V[0::2][:k + 1]

        # Halve n and iterate until n == 0
        n = n >> 1

    # The n-th coefficient of the original P(t) / Q(t) is now equal to the 0-th
    # coefficient of P(t) / Q(t) after log_2(n) iterations.
    return P[0] // Q[0]

