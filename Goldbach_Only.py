import math
import numpy as np

# 24 classes coprime to 90
ALL_CLASSES = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89]

def epoch_size(h):
    return 90*h*h - 12*h + 1

# Operator derivation exactly as in your papers
def get_operators(k):
    operators = []
    seen = set()
    for z in ALL_CLASSES:
        try:
            o = (k * pow(z, -1, 90)) % 90
            if o not in ALL_CLASSES: continue
            pair = tuple(sorted([z,o]))
            if pair in seen: continue
            seen.add(pair)
            z_eff = 91 if z==1 else z
            o_eff = 91 if o==1 else o
            l = 180 - (z_eff + o_eff)
            m = 90 - (z_eff + o_eff) + (z_eff*o_eff - k)//90
            operators.append((l,m,z_eff))
            if z != o:
                operators.append((l,m,o_eff))
        except: continue
    return operators

def hole_primes_up_to_h(h):
    N = epoch_size(h)
    all_primes = set()
    xmax = int(math.sqrt(250*N/90)) + 10
    for k in ALL_CLASSES:
        amplitude = np.zeros(N, dtype=bool)
        ops = get_operators(k)
        for l,m,z in ops:
            for x in range(1,xmax+1):
                y = 90*x*x - l*x + m
                if y >= N: continue
                p = z + 90*(x-1)
                if p <= 0: continue
                idx = int(y)
                while idx < N:
                    amplitude[idx] = True
                    idx += p
        holes_n = np.where(~amplitude)[0]
        primes = {90*n + k for n in holes_n if 90*n + k > 5}
        all_primes.update(primes)
    return sorted(all_primes)

def goldbach_multiplicity_in_range(E_start, E_end, step, h_values=[50,100,150,200]):
    print("Goldbach Variant (Helkenberg–Sieve) – Multiplicity Growth")
    print("OEIS A002375 – Number of ways 2n = p+q (p≤q odd primes): https://oeis.org/A002375\n")
    results = {}
    for h in h_values:
        print(f"--- Epoch h={h} (N≈{epoch_size(h):,}) ---")
        primes = hole_primes_up_to_h(h)
        print(f"Total hole primes extracted: {len(primes)}")
        mult = []
        for E in range(E_start, E_end+1, step):
            count = 0
            seen = set()
            for p in primes:
                if p >= E//2 + 1: break
                q = E - p
                if q in primes:
                    pair = tuple(sorted([p,q]))
                    if pair not in seen:
                        seen.add(pair)
                        count += 1
            mult.append(count)
        avg = np.mean(mult)
        minimum = min(mult)
        print(f"Even range {E_start}–{E_end}: average {avg:.2f} ways, minimum {minimum}")
        results[h] = (avg, minimum)
    return results

# Example run (adjust range/h as desired)
if __name__ == "__main__":
    goldbach_multiplicity_in_range(10000, 20000, 100, h_values=[50,80,110,140])