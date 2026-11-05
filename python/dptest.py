from mpmath import nint, fabs, floor, log

dist = lambda x: fabs(x - nint(x))

def DP_test(M, gamma, mu, log_A, log_B):
    target = 6 * M
    
    eta = gamma
    q_old = 1
    q = 0

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

    N = (log_A + log(q) - log(epsilon)) / log_B

    return {"q": q, "epsilon": epsilon, "N": N}