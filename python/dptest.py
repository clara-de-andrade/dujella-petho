#!/usr/bin/python3

"""
This modules contains a function dptest implementing the Dujella-Petho
reduction. Given an integer M > 0 and irrational numbers gamma and mu,
it computes the least denominator q of a convergent p/q of the continued
fraction of gamma such that two conditions are met:

    1) q > 6 * M,
    2) epsilon > 0;

where epsilon = ||q * mu|| - M * ||q * gamma||, with ||x|| denoting the
distance of x from the nearest integer to x. Satisfied these conditions, the
Dujella-Petho reduction lemma states there are no positive integers m, n, k
solving the inequality

    0 < |m * gamma - k + mu| < A * B**n

with m <= M and n >= N, where N = log(A * q / epsilon) / log(B). The function
DP_test also computes N and returns a dictionary containing q, epsilon and N.
"""

from mpmath import nint, fabs, floor, log
from mpmath import mpf


dist = lambda x: fabs(x - nint(x))

def dptest(M: int, gamma: mpf, mu: mpf, log_A: mpf, log_B: mpf):
    """
    Computes q, epsilon, verifies the conditions for the Dujella-Petho
    reduction lemma, and returns q, epsilon and the reduced upper limit.
    """

    target = 6 * M
    
    # Initial values of eta, q_old and q
    eta = gamma
    q_old = 1
    q = 0

    # Computes q = q_i by the reccurence relation
    #
    #    q_i = a_i * q_(i - 1) + q_(i - 2),
    #
    # with q_(-2) = 1 and q_(-1) = 0, where a_i is the i-th coefficient of
    # the continued fraction for gamma, i.e.
    #
    #   gamma = 1 / (a_0 + 1 / (a_1 + 1 / (a_2 + 1 / ...))).
    #
    # In particular, a_i is computed iteratively by putting
    #
    #   eta_0 = gamma, eta_(i + 1) = 1 / (eta_i - a_i), a_i = floor(eta_i).

    while True:
        a = floor(eta)
        eta = 1 / (eta - a)
        q_old, q = q, q * a + q_old
        if not q > target:
            continue

        epsilon = dist(q * mu) - M * dist(q * gamma)
        if not epsilon > 0:
            continue

        break

    # Once the conditions are satisfied, compute the reduced upper limit
    N = (log_A + log(q) - log(epsilon)) / log_B

    return {"q": q, "epsilon": epsilon, "N": N}

