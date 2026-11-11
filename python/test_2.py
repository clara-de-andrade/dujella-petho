from mpmath import mp, mpf, nstr, log, floor, sqrt
from math import floor as m_floor
from math import ceil as m_ceil

from dptest import dptest

import pickle


l_min = 2
k_min = 2
q_min = 3

a_min = 1
a_max = 9

def prog_bar(i, i_min, i_max):
    prog = int(100 * (i - i_min + 1) / (i_max - i_min + 1))
    return f"[{"#" * m_floor(prog / 2)}{" " * m_ceil((100 - prog) / 2)}] {prog}%   "


mp.prec = 1 << 12
q_max = 10

if __name__ == "__main__":
    print("Calculando cotas iniciais...")

    k_max = {}
    for q in range(q_min, q_max + 1):
        k_max[q] = floor(mpf("1.6e19") * (log(q + 1) ** 4))

    l_max = {}
    for q in range(q_min, q_max + 1):
        l_max[q] = floor(mpf("1.2e95") * (log(q + 1) ** 21))
        # l_max[q] = floor(mpf("3.84e34") * (q ** 5) * (log(q + 1) ** 21))

    print("PRONTO!\n")

    print("Calculando alpha e Gamma...")
    
    alpha = {}
    Gamma = {}
    for q in range(q_min, q_max + 1):
        sqrtD = sqrt((q - 1) ** 2 + 4)

        alpha[q] = (q + 1 + sqrtD) / 2
        Gamma[q] = 1 + (q - 1) / sqrtD

        print(prog_bar(q, q_min, q_max), end="\r")
    print("\nPRONTO!\n")
            
    print("Calculando epsilon, q e N...")

    tests = []
    for q in range(q_min, q_max + 1):
        gamma = log(10) / log(alpha[q])
        log_A = log(2 * q)
        log_B = log(q) / 2
        for a in range(a_min, a_max + 1):
            mu = (log(a / (9 * Gamma[q]))) / log(alpha[q])

            tests.append(dptest(l_max[q], gamma, mu, log_A, log_B))
            tests[-1].update({"(q, a)": (q, a)})

        print(prog_bar(q, q_min, q_max), end="\r")
    print("\nPRONTO!")

    with open("test_2.pkl", "wb") as f:
        pickle.dump(tests, f)