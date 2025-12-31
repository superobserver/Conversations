import math
import cmath
from sympy import lcm, isprime

def drLD(x, l, m, z, listvar):
    y = 90 * (x * x) - l * x + m
    if y >= len(listvar) or y < 0:
        return
    listvar[int(y)] += 1
    p = z + (90 * (x - 1))
    if p == 0:
        return
    for n in range(1, int(((len(listvar) - 1 - y) / p) + 1)):
        next_y = y + (p * n)
        if next_y >= len(listvar):
            break
        listvar[int(next_y)] += 1

def harmonic_disharmonic_sieve(h=10, residue_k=17):
    epoch = 90 * (h * h) - 12 * h + 1
    a = 90
    b = -300
    c = 250 - epoch
    d = (b ** 2) - (4 * a * c)
    sol2 = (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = sol2.real

    list_amp = [0] * int(epoch + 100)

    operators = [
        (72, -1, 17), (72, -1, 91),
        (108, 29, 19), (108, 29, 53),
        (72, 11, 37), (72, 11, 71),
        (18, 0, 73), (18, 0, 89),
        (102, 20, 11), (102, 20, 67),
        (138, 52, 13), (138, 52, 29),
        (102, 28, 31), (102, 28, 47),
        (48, 3, 49), (48, 3, 83),
        (78, 8, 23), (78, 8, 79),
        (132, 45, 7), (132, 45, 41),
        (78, 16, 43), (78, 16, 59),
        (42, 4, 61), (42, 4, 77)
    ]

    # Harmonic marking (original)
    periods = set()
    for x in range(1, int(new_limit) + 1):
        for l, m, z in operators:
            drLD(x, l, m, z, list_amp)
            p = z + 90 * (x - 1)
            if p > 1:
                periods.add(p)

    list_amp = list_amp[:epoch]

    # LCM matrix (24x24)
    lcm_matrix = [[1] * 24 for _ in range(24)]
    op_periods = [[] for _ in operators]
    for op_idx, (l, m, z) in enumerate(operators):
        for x in range(1, int(new_limit) + 1):
            p = z + 90 * (x - 1)
            if p > 1:
                op_periods[op_idx].append(p)
    for i in range(24):
        for j in range(i, 24):
            if op_periods[i] and op_periods[j]:
                lcm_val = lcm(op_periods[i][0], op_periods[j][0])  # Proxy first for small h
            else:
                lcm_val = 1
            lcm_matrix[i][j] = lcm_matrix[j][i] = lcm_val

    # Q = LCM of all unique periods
    Q = 1
    for p in periods:
        Q = lcm(Q, p)

    # Anti-chain generator: y coprime to Q (disharmonics)
    anti_chain_y = []
    for y in range(epoch):
        n = 90 * y + residue_k
        if math.gcd(n, Q) == 1 and list_amp[y] == 0:  # Verify with amplitude
            anti_chain_y.append(y)

    # Output for comparison
    harmonic_holes = [y for y, amp in enumerate(list_amp) if amp == 0]
    print("Harmonic holes (primes):", len(harmonic_holes))
    print("Anti-chain holes (predicted primes):", len(anti_chain_y))
    print("Match:", set(anti_chain_y) == set(harmonic_holes))
    print("Sample LCM matrix row 0:", lcm_matrix[0][:5])
    return anti_chain_y, lcm_matrix

if __name__ == "__main__":
    harmonic_disharmonic_sieve(h=10)
