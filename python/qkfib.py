#!/usr/bin/python3

"""
This module contains a function qkfib, which computes the n-th term of the
(q, k)-generalised Fibonacci sequence, for brevity denoted simply by F_n, where

    F_n = q * F_(n - 1) + F_(n - 2) + ... + F_(n - k),                 (n >= 2)
    F_(-k + 2) = ... = F_0 = 0,
    F_1 = 1.

The implementation of qkfib uses the generating function of F_n, given by

    f(q, k)(t) = t / (1 - qt - t**2 - ... - t**k)
               = sum(F(q, k, i) * t**i for i >= 0)

for all t in its radius of convergence. (In what follows, q >= 1 and k >= 2.)
"""

from math import log2, floor


def polymul(A_in, B_in):
    """
    Computes the product of two polynomials A(t) = a_0 + ... + a_m * t**m and
    B = b_0 + ... + b_n * t**n, where m, n >= 0, using Karatsuba's algorithm.
    """

    # Initialize A and B
    A = list(A_in)
    B = list(B_in)

    # Force deg(A) >= deg(B)
    if len(A) < len(B):
        A, B = B, A

    # Set d to smallest power of 2 larger than deg A
    d = len(A)
    if (d & (d - 1)) != 0:
        d = 1 << (floor(log2(d)) + 1)

    # Force deg(A) = deg(B) = d by padding A and B with zeroes
    A = A + [0] * (d - len(A))
    B = B + [0] * (d - len(B))

    # Trivial case
    if d == 1:
        return [A[0] * B[0]]

    # m = d / 2
    m = d >> 1
    # A = A0 + (t**m) * A1
    A0, A1 = A[:m], A[m:]
    # B = B0 + (t**m) * B1
    B0, B1 = B[:m], B[m:]

    # By the above equations for A and B, we have that
    #
    #    A * B = (A0 + (t**m) * A1) * (B0 + (t**m) * B1)
    #          = A0 * B0 + (t**m) * (A1 * B0 + A0 * B1) + (t**d) (A1 * B1).
    #
    # In particular, we also have that
    #
    #    z3 = A1 * B0 + A0 * B1
    #       = (A0 + A1) * (B0 + B1) - A0 * B0 - A1 * B1.
    #
    # Therefore, we have reduced the problem of computing one product of
    # polynomials of degree d to that of computing just three products of
    # polynomials of degree d / 2.

    # z0 = A0 * B0
    z0 = polymul(A0, B0)
    # z1 = (A0 + A1) * (B0 + B1)
    z1 = polymul(
        [A0[i] + A1[i] for i in range(m)],
        [B0[i] + B1[i] for i in range(m)]
    )
    # z2 = A1 * B1
    z2 = polymul(A1, B1)

    # z3 = A1 * B0 + A0 * B1
    #    = (A0 + A1) * (B0 + B1) - (A0 * B0 + A1 * B1).
    #    = z1 - z0 - z2
    z3 = [z1[i] - (z0[i] + z2[i]) for i in range(len(z1))]

    # Initialize M, with deg M = 2 * d, and set to [0, ..., 0].
    M = [0] * (d << 1)
    
    # Setting M to z0 + (t**m) * z3 + (t**d) * z2 = A * B

    # M += z0
    for i in range(len(z0)):
        M[i] += z0[i]
    # M += (t**m) * z3
    for i in range(len(z3)):
        M[m + i] += z3[i]
    # M += (t**d) * z2
    for i in range(len(z2)):
        M[d + i] += z2[i]

    # Remove leading zero coefficients of M
    M = M[:len(A_in) + len(B_in) - 1]

    return M


def qkfib(q, k, n):
    """
    Computes the n-th term of the (q, k)-generalised Fibonacci sequence, as the
    n-th coefficient of the formal power series for its generating function,
    using the Bostan-Mori algorithm. (Algorithm 1 of arXiv:2008.08822v1.)
    """

    # F(q, k, n) = 0 for all -k + 2 <= n <= 0
    if n < 1:
        return 0
    
    # P(t) = t
    P = [0, 1] + [0] * (k - 2)
    # Q(t) = 1 - qt - t**2 - ... - t**k
    Q = [1,-q] + [-1] * (k - 1)

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
        #     m = floor(n / 2)
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

