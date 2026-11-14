from mpmath import mp, mpf, nstr, log, floor
from math import floor as m_floor
from math import ceil as m_ceil

from polynom import Poly

from dptest import dptest
from findroot import findroot

import pickle


m_min = 3
k_min = 3

def prog_bar(i, i_min, i_max):
    prog = int(100 * (i - i_min + 1) / (i_max - i_min + 1))
    return f"[{"#" * m_floor(prog / 2)}{" " * m_ceil((100 - prog) / 2)}] {prog}%   "


mp.prec = 4096
m_max = 1394

if __name__ == "__main__":
    print("Calculando cotas iniciais...")

    s_max = floor(mpf("4.7e14") * (m_max ** 5) * (log(m_max) ** 4))  
    print(f"s <= {nstr(s_max, 3)}")

    n_max = m_max * s_max + 2
    print(f"n <= {nstr(n_max, 3)}")

    k_max = int(floor(log(s_max)))
    print(f"k <= {k_max}")
    print("PRONTO!\n")
    cache = {"k_max": k_max, "n_max": n_max, "s_max": s_max, "m_max": m_max}

    print("Calculando F...")

    F = {}
    for k in range(k_min, k_max + 1):
        F[k, 0] = 0
        F[k, 1] = 1
        prev = [0] * (k - 1) + [1]
        for m in range(2, m_max + 2):
            F[k, m] = sum(prev)
            prev = prev[1:] + [F[k, m]]

        print(prog_bar(k, k_min, k_max), end="\r")
    print("\nPRONTO!\n")
    cache.update({"F": F})

    print("Calculando alpha e g...")

    A = mpf("2.01"))
    B = mpf("1.65") ** (1 / mpf(m_max))
    cache.update({"A": mpf("2.01"), "B": mpf("1.65") ** (1 / mpf(m_max))})

    alpha = {}
    log_alpha = {}
    g = {}
    log_g = {}
    for k in range(k_min, k_max + 1):
        P = Poly(k, [-1] * k + [1])

        alpha[k] = findroot(P, 2)
        g[k] = (alpha[k] - 1) / (2 + (k + 1) * (alpha[k] - 2))

        log_alpha[k] = log(alpha[k])
        log_g[k] = log(g[k])

        print(prog_bar(k, k_min, k_max), end="\r")
    print("\nPRONTO!\n")
    cache.update({"alpha": alpha, "g": g})
            
    print("Calculando epsilon, q e N...")

    gamma = {}
    mu = {}
    q = {}
    epsilon = {}
    N = {}
    for m in range(m_min, m_max + 1):
        for k in range(k_min, min(m, k_max) + 1):
            log_F = log(F[k, m + 1])

            gamma[k, m] = log_alpha[k] / log_F
            mu[k, m] = (log_alpha[k] - log_g[k]) / log_F
            
            test = dptest(n_max, gamma[k, m], mu[k, m], A, B)
            N[k, m] = test["N"]
            q[k, m] = test["q"]
            epsilon[k, m] = test["epsilon"]

        print(prog_bar(m, m_min, m_max), end="\r")
    print("\nPRONTO!")
    cache.update({"gamma": gamma, "mu": mu, "q": q, "epsilon": epsilon, "N": N})

    print("Calculando novas cotas...")
    n_max = floor(max(cache["N"].values()))
    print(f"n <= {nstr(n_max, 3)}")

    s_max = floor((n_max - 2) / (m_min - 1))
    print(f"s <= {nstr(s_max, 3)}")
    
    k_max = int(floor(log(s_max)))
    print(f"k <= {k_max}")
    print("PRONTO!\n")

    cache.update({"n_max_new": n_max, "s_max_new": s_max, "k_max_new": k_max})

    gamma = {}
    mu = {}
    q = {}
    epsilon = {}
    N = {}
    for m in range(m_min, m_max + 1):
        for k in range(k_min, min(m, k_max) + 1):
            log_F = log(F[k, m + 1])

            gamma[k, m] = log_alpha[k] / log_F
            mu[k, m] = (log_alpha[k] - log_g[k]) / log_F
            
            test = dptest(n_max, gamma[k, m], mu[k, m], log_A, log_B)
            q[k, m] = test["q"]
            epsilon[k, m] = test["epsilon"]
            N[k, m] = test["N"]

        print(prog_bar(m, m_min, m_max), end="\r")
    print("\nPRONTO!")
    cache.update({"gamma": gamma, "mu": mu, "q": q, "epsilon": epsilon, "N": N})

    with open("test_1.pkl", "wb") as f:
        pickle.dump(cache, f)