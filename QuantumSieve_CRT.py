import argparse
import math
from sympy import isprime

def mod_inverse(a, m):
    try:
        return pow(a, -1, m)
    except ValueError:
        return None

def get_residues():
    return [r for r in range(1,90) if math.gcd(r,90)==1]  # Auto-generate 24 primitives (digital root/last digit rules implicit)

def get_pairs(k, residues):
    pairs = []
    seen = set()
    for r in residues:
        inv = mod_inverse(r, 90)
        if inv is not None:
            s = (k * inv) % 90
            if s in residues and s >= r and (r, s) not in seen:
                seen.add((r, s))
                pairs.append((r, s))
    return pairs  # 12 pairs

def adjust_reps(r, s, k):
    r_adj, s_adj = r, s
    n0 = (r_adj * s_adj - k) // 90
    num = 90 * n0 + k
    if n0 <= 0 or (num > 1 and isprime(num)):
        # Adjust the one closer to 1 or add 90 to smaller
        if r_adj < s_adj:
            r_adj += 90
        else:
            s_adj += 90
        n0 = (r_adj * s_adj - k) // 90
    return r_adj, s_adj

def get_coeffs(r, s, k):
    r_adj, s_adj = adjust_reps(r, s, k)
    l = 180 - (r_adj + s_adj)
    n0 = (r_adj * s_adj - k) // 90
    m = n0 + 90 - l
    return l, m

def drLD(x, l, m, z, listvar, limit):
    y = 90 * x**2 - l * x + m
    if y >= limit or y < 0:
        return
    listvar[y] += 1
    p = z + 90 * (x - 1)
    for step in range(1, (limit - y) // p + 1):
        next_pos = y + p * step
        if next_pos < limit:
            listvar[next_pos] += 1

def run_sieve(h, k):
    residues = get_residues()
    pairs = get_pairs(k, residues)
    epoch = 90 * h**2 - 12 * h + 1
    limit = epoch
    listvar = [0] * (limit + 100)
    new_limit = int(math.sqrt(limit / 90)) + 5  # Approx x iters

    for x in range(1, new_limit + 1):
        for r, s in pairs:
            l, m = get_coeffs(r, s, k)
            drLD(x, l, m, r, listvar, limit)
            drLD(x, l, m, s, listvar, limit)

    listvar = listvar[:limit]
    marks = sum(listvar)
    primes = sum(1 for amp in listvar if amp == 0)
    x_iters = new_limit
    operators = 24 * x_iters  # 12 pairs * 2

    print(f"Epoch: {epoch}")
    print(f"Marks: {marks}")
    print(f"x iters: {x_iters}")
    print(f"Operators: {operators}")
    print(f"Quantity Primes: {primes}")
    return primes

def run_twins(h, k1, k2):
    primes1 = run_sieve(h, k1)  # Run for k1
    primes2 = run_sieve(h, k2)  # Run for k2
    # To find common holes (twins), intersect unmarked n, but for reproduction, print separate

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Deterministic Quadratic Sieve')
    parser.add_argument('--k', type=int, required=True, help='Residue class k')
    parser.add_argument('--h', type=int, required=True, help='Epoch parameter h')
    parser.add_argument('--k_twin', type=int, default=None, help='Twin class k2 for twins')
    args = parser.parse_args()

    if args.k_twin:
        run_twins(args.h, args.k, args.k_twin)
    else:
        run_sieve(args.h, args.k)