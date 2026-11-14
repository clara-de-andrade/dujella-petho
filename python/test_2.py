from mpmath import mp, mpf, nstr, log, floor, sqrt
from math import floor as m_floor
from math import ceil as m_ceil

from dptest import dptest

import pickle


l_min = 2
k_min = 2
q_min = 2

a_min = 1
a_max = 9

def prog_bar(i, i_min, i_max):
    prog = int(100 * (i - i_min + 1) / (i_max - i_min + 1))
    return f"[{"#" * m_floor(prog / 2)}{" " * m_ceil((100 - prog) / 2)}] {prog}%   "


mp.prec = 1 << 12
q_max = 1000

if __name__ == "__main__":
    print("Calculando cotas iniciais...")

    l_max = {}
    for q in range(q_min, q_max + 1):
        l_max[q] = floor(mpf("1.2e95") * (log(q + 1) ** 21))
        print(prog_bar(q, q_min, q_max), end="\r")
    print("\nPRONTO!\n")
    cache = {"l_max": l_max, "q_max": q_max}

    print("Calculando alpha e Gamma...")
    
    alpha = {}
    Gamma = {}
    for q in range(q_min, q_max + 1):
        sqrtD = sqrt((q - 1) ** 2 + 4)

        alpha[q] = (q + 1 + sqrtD) / 2
        Gamma[q] = 1 + (q - 1) / sqrtD

        print(prog_bar(q, q_min, q_max), end="\r")
    print("\nPRONTO!\n")
    cache = {"alpha": alpha, "Gamma": Gamma}
            
    print("Calculando epsilon, q e N...")

    A = {q: 2 * q for q in range(q_min, q_max + 1)}
    B = {q: mpf(q) ** 0.5 for q in range(q_min, q_max + 1)}
    cache.update({"A": A, "B": B})

    gamma = {}
    mu = {}
    q_conv = {}
    epsilon = {}
    N = {}
    for q in range(q_min, q_max + 1):
        gamma[q] = log(10) / log(alpha[q])
        for a in range(a_min, a_max + 1):
            mu[q, a] = (log(a / (9 * Gamma[q]))) / log(alpha[q])

            test = dptest(l_max[q], gamma[q], mu[q, a], A[q], B[q])
            q_conv[q, a] = test["q"]
            epsilon[q, a] = test["epsilon"]
            N[q, a] = test["N"]

        print(prog_bar(q, q_min, q_max), end="\r")
    print("\nPRONTO!")
    cache.update({"gamma": gamma, "mu": mu, "q": q_conv, "epsilon": epsilon, "N": N})

    with open("test_2.pkl", "wb") as f:
        pickle.dump(cache, f)