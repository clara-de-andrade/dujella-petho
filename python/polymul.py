from math import log2, floor
    
def polymul(A_in, B_in):
    A = list(A_in)
    B = list(B_in)

    if len(A) < len(B):
        A, B = B, A

    size = len(A)
    if (size & (size - 1)) != 0:
        size = 1 << (floor(log2(size)) + 1)

    A = A + [0] * (size - len(A))
    B = B + [0] * (size - len(B))

    if size == 1:
        return [A[0] * B[0]]
    half = size >> 1

    A0, A1 = A[:half], A[half:]
    B0, B1 = B[:half], B[half:]

    z0 = polymul(A0, B0)
    z1 = polymul(
        [A0[i] + A1[i] for i in range(half)],
        [B0[i] + B1[i] for i in range(half)]
    )
    z2 = polymul(A1, B1)
    z3 = [z1[i] - z0[i] - z2[i] for i in range(len(z1))]

    M = [0] * (size << 1)
    for i in range(len(z0)):
        M[i] += z0[i]

    for i in range(len(z3)):
        M[half + i] += z3[i]

    for i in range(len(z2)):
        M[size + i] += z2[i]

    M = M[:len(A_in) + len(B_in) - 1]

    return M