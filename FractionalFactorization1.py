import math
import numpy as np
import matplotlib.pyplot as plt
from fractions import Fraction

# 24 residue classes from the documents
ALL_CLASSES = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89]

# Function to get operators (from your sieve logic)
def get_operators(k):
    operators = []
    seen = set()
    for z in ALL_CLASSES:
        try:
            o = (k * pow(z, -1, 90)) % 90
            if o not in ALL_CLASSES: continue
            pair = tuple(sorted([z, o]))
            if pair in seen: continue
            seen.add(pair)
            z_eff = 91 if z == 1 else z
            o_eff = 91 if o == 1 else o
            l = 180 - (z_eff + o_eff)
            m = 90 - (z_eff + o_eff) + (z_eff * o_eff - k) // 90
            operators.append((l, m, z_eff, z))
            if z != o:
                operators.append((l, m, o_eff, o))
        except ValueError:
            continue
    return operators

# y_x and p_x
def y_x(x, l, m):
    return 90 * x**2 - l * x + m

def p_x(x, z_eff):
    return z_eff + 90 * (x - 1)

# Max x where y_x < n
def max_x_for_n(n):
    return int(math.sqrt(n / 90)) + 1

# Get fractional parts for failed tests
def get_fractional_parts(n, k):
    ops = get_operators(k)
    fractions = []
    for l, m, z_eff, z in ops:
        max_x = max_x_for_n(n)
        for x in range(1, max_x + 1):
            y = y_x(x, l, m)
            if y >= n: continue
            p = p_x(x, z_eff)
            if p <= 0: continue
            raw_b = (n - y) / p
            frac = raw_b - math.floor(raw_b)
            if frac != 0:  # Failed (non-integer b)
                fractions.append((x, frac))
    return fractions

# Check if n is hole (no chaining) or chained
def is_hole(n, k):
    ops = get_operators(k)
    max_x = max_x_for_n(n)
    for l, m, z_eff, z in ops:
        for x in range(1, max_x + 1):
            y = y_x(x, l, m)
            if y >= n: continue
            p = p_x(x, z_eff)
            if p <= 0: continue
            raw_b = (n - y) / p
            if raw_b == math.floor(raw_b) and raw_b >= 0:
                return False  # Chained
    return True  # Hole

# Main simulation
n_start = 10000
n_end = 10100
hole_tuples = []
chained_tuples = []

for n in range(n_start, n_end + 1):
    k = 11  # Example class; loop over ALL_CLASSES for full 24
    fractions = get_fractional_parts(n, k)
    if is_hole(n, k):
        hole_tuples.append((n, fractions))
    else:
        chained_tuples.append((n, fractions))

# Sample output
print("Sample Hole Tuples (n, [(x, frac), ...]):")
for tup in hole_tuples[:3]:
    print(tup)

print("\nSample Chained Tuples:")
for tup in chained_tuples[:3]:
    print(tup)

# Graph fractional parts vs x
fig, ax = plt.subplots(2, 1, figsize=(12, 8))

# Holes
for n, fracs in hole_tuples:
    xs = [t[0] for t in fracs]
    fs = [t[1] for t in fracs]
    ax[0].plot(xs, fs, marker='o', alpha=0.6, label=f'n={n}')
ax[0].set_title('Fractional Parts for Holes (Primes/Holes)')
ax[0].set_xlabel('x')
ax[0].set_ylabel('Fractional Part')
ax[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Chained
for n, fracs in chained_tuples:
    xs = [t[0] for t in fracs]
    fs = [t[1] for t in fracs]
    ax[1].plot(xs, fs, marker='o', alpha=0.6, label=f'n={n}')
ax[1].set_title('Fractional Parts for Chained (Composites)')
ax[1].set_xlabel('x')
ax[1].set_ylabel('Fractional Part')
ax[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.show()  # Or save to file for analysis