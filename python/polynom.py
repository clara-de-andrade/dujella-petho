from dataclasses import dataclass
from math import log2, floor
from mpmath import mpf

@dataclass
class Poly:
    deg: int = 0
    coeffs: list = [0]

    def __add__(self, other):
        # TODO
        pass

    def __mul__(self, other):
        # TODO
        pass

    def diff():
        # TODO
        pass

    def eval_at(x):
        # TODO
        pass


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