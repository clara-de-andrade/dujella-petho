#!/usr/bin/python3

# TODO: Document module
"""
"""

from polynom import Poly


def nthterm(P: Poly, Q: Poly, n : int) -> Number:
    """
    ...
    using the Bostan-Mori algorithm. (Algorithm 1 of arXiv:2008.08822v1.)
    """

    if not P.deg < Q.deg:
        raise ValueError

    while n > 0:
        # Q_minus(t) = Q(-t)
        Q_minus = Poly(k, [Q[i] * (-1) ** (i % 2) for i in range(k + 1)])
        # U(t) = P(t) * Q(-t)
        U = P * Q_minus
        # V(t) = Q(t) * Q(-t)
        V = Q * Q_minus

        # Note that, since P(t) / Q(t) = U(t) / V(t**2), if we let
        #
        #     U(t) = A(t**2) + t * B(t**2),
        #     m = floor(k / 2)
        #
        # then the n-th coefficient of P(t) / Q(t) is equal to the m-th
        # coefficient of A(t)/V(t) if n is even, and equal to the m-th
        # coefficient of B(t)/V(t) if n is odd.

        if n & 1 == 0:
            P = U[0::2][:k]
        else:
            P = U[1::2][:k]
        Q = V[0::2][:k + 1]

        # Halve n and iterate until n == 0
        n = n >> 1

    # The n-th coefficient of the original P(t) / Q(t) is now equal to the 0-th
    # coefficient of P(t) / Q(t) after log_2(n) iterations.
    return P[0] // Q[0]

def genfunc(k: int, coeffs: list, inits: list, n: int):
    """
    ...
    """

    if n >= 0 and n < k:
        return inits[n]

    P = Poly(
        k - 1,
        [inits[n] - sum(
            coeffs[i - 1] * inits[n - i]
            for i in range(1, n + 1))
         for n in range(k)]
    )
    Q = -Poly(k, [-1] + coeffs)

    return P, Q
