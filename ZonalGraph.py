import math
import cmath
from math import gcd
from functools import reduce
import matplotlib.pyplot as plt

def lcm(a, b):
    return abs(a * b) // gcd(a, b) if a and b else 0

def lcm_list(nums):
    return reduce(lcm, nums, 1)

def get_operators(k):
    R = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89]
    operators = []
    seen_pairs = set()
    for z in R:
        try:
            o = (k * pow(z, -1, 90)) % 90
            if o not in R:
                continue
            pair = tuple(sorted([z, o]))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            z_eff = 91 if z == 1 else z
            o_eff = 91 if o == 1 else o
            l = 180 - (z_eff + o_eff)
            m = 90 - (z_eff + o_eff) + (z_eff * o_eff - k) // 90
            operators.append((l, m, z_eff, k))
            if z != o:
                operators.append((l, m, o_eff, k))
        except ValueError:
            continue
    return operators

def zonal_sieve_segment(h, classes, segment_start, segment_end):
    a = 90
    b = -300
    c = 250 - (90 * (h * h) - 12 * h + 1) 
    d = (b**2) - (4*a*c)
    sol2 = (-b + cmath.sqrt(d)) / (2*a)
    new_limit = int(sol2.real) + 1
    segment_width = segment_end - segment_start
    list_amp = [0] * segment_width 
    all_ops = []
    for k in classes:
        all_ops.extend(get_operators(k))
    for x in range(1, new_limit + 1):
        for l_val, m_val, z, primitive in all_ops:
            y0 = 90 * x * x - l_val * x + m_val
            p = z + 90 * (x - 1)
            if p <= 0 or y0 >= segment_end:
                continue
            if y0 < segment_start:
                diff = segment_start - y0
                n = math.ceil(diff / p)
                current = y0 + n * p
            else:
                current = y0
            while current < segment_end:
                if current >= segment_start:
                    idx = int(current - segment_start)
                    list_amp[idx] += 1
                current += p
    return list_amp 

R = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89]
segment_start = 0
segment_end = 5000
h = int(math.sqrt(segment_end * 4)) + 10 

cumulatives = {}
for r in R:
    classes = [r]
    amp_segment = zonal_sieve_segment(h, classes, segment_start, segment_end)
    cum = [0]
    for a in amp_segment:
        cum.append(cum[-1] + (1 if a == 0 else 0))
    cumulatives[r] = cum
    print(f"For {r}: holes {cum[-1]}")

fig, ax = plt.subplots(figsize=(12,8))
for r, cum in cumulatives.items():
    ax.plot(range(len(cum)-1), cum[1:], label=str(r))
ax.legend(title="Residue class")
ax.set_xlabel("k (up to 4999)")
ax.set_ylabel("Cumulative primes in 90k + r")
ax.set_title("Race for highest prime count in residue classes mod 90")
plt.savefig('race.png')
print("Plot saved to race.png")




#Real life Chun-Li. #kickboxing
