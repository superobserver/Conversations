#!/usr/bin/env python3
import cmath
import math
from functools import reduce

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    return a * b // gcd(a, b)

def lcmm(*args):
    return reduce(lcm, args)

h = 10
epoch = 90 * (h * h) - 12 * h + 1
limit = epoch

a, b, c = 90, -300, 250 - limit
new_limit = (-b + cmath.sqrt(b**2 - 4*a*c)) / (2*a)

is_marked = [False] * int(limit + 100)

params = [
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

for x in range(1, int(new_limit.real)):
    ys = []
    ps = []
    for l, m, z in params:
        y = 90 * (x * x) - l * x + m
        p = z + 90 * (x - 1)
        ys.append(int(y))
        ps.append(p)
    
    lcm_p = lcmm(*ps)
    
    for i in range(24):
        y = ys[i]
        p = ps[i]
        if y >= limit: continue
        is_marked[y] = True
        pos = y + p
        while pos < limit:
            is_marked[pos] = True
            pos += p
    
    # LCM optimization: mark multiples of LCM starting from min y
    min_y = min(ys)
    pos = min_y + lcm_p
    while pos < limit:
        is_marked[pos] = True  # Since all operators align here
        pos += lcm_p

is_marked = is_marked[:limit]

holes = [i for i, marked in enumerate(is_marked) if not marked]

primes_small = [90 * n + 11 for n in holes if 90 * n + 11 > 1]
primes_large = [90 * n + 13 for n in holes if 90 * n + 13 > 1]

print(f"Holes: {len(holes)}")
print("First 10 holes:", holes[:10])
print("First 10 smaller primes:", primes_small[:10])
print("First 10 larger primes:", primes_large[:10])