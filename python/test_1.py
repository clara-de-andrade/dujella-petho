from mpmath import mp, mpf, nstr, log, floor
from math import floor as m_floor
from math import ceil as m_ceil

from polynom import Poly

from dptest import dptest
from nthterm import nthterm
from findroot import findroot

import pickle


m_min = 3
k_min = 3

def prog_bar(i, i_min, i_max):
    prog = int(100 * (i - i_min + 1) / (i_max - i_min + 1))
    return f"[{"#" * m_floor(prog / 2)}{" " * m_ceil((100 - prog) / 2)}] {prog}%   "

def F(k, n):
    P = Poly(1, [0, 1])
    Q = Poly(k, [1] + [-1] * k)
    return nthterm(P, Q, n)


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

    print("Calculando alpha e g...")
    
    alpha = {}
    g = {}
    for k in range(k_min, k_max + 1):
        P = Poly(k, [-1] * k + [1])

        alpha[k] = findroot(P, 2)
        g[k] = (alpha[k] - 1) / (2 + (k + 1) * (alpha[k] - 2))

        print(prog_bar(k, k_min, k_max), end="\r")
    print("\nPRONTO!\n")
            
    print("Calculando epsilon, q e N...")

    log_A = log(mpf("2.01"))
    log_B = log(mpf("1.65"))

    tests = []
    for m in range(m_min, m_max + 1):
        for k in range(k_min, min(m, k_max) + 1):
            gamma = log(alpha[k]) / log(F(k, m + 1))
            mu = log(alpha[k] / g[k]) / log(F(k, m + 1))
            
            tests.append(dptest(n_max, gamma, mu, log_A, log_B))
            tests[-1]["N"] *= m_max
            tests[-1].update({"(m, k)": (m, k)})

        print(prog_bar(m, m_min, m_max), end="\r")
    print("\nPRONTO!")

    with open("test_1.pkl", "wb") as f:
        pickle.dump(tests, f)