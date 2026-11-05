from polymul import polymul

def F(q, k, n):
    if n == 0:
        return 0
    if k == 1:
        return q ** (n - 1)
    
    P = [0, 1] + [0] * (k - 2)
    Q = [1,-q] + [-1] * (k - 1)

    while n > 0:
        Q_minus = [Q[i] if i % 2 == 0 else -Q[i] for i in range(k + 1)]
        U = polymul(P, Q_minus)

        if n % 2 == 0:
            P = U[0::2][:k]
        else:
            P = U[1::2][:k]

        V = polymul(Q, Q_minus)
        Q = V[0::2][:k + 1]

        n = n >> 1

    return P[0] // Q[0]

