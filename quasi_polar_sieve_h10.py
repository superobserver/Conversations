#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quasi-polar visualisation of the quadratic twin-prime sieve (h = 10)

Re-implements NewGrok17.py exactly, then draws the markings on a polar grid:
    • 24 rays  →  one per drLD call
    • radius  →  address i  (scaled so the outermost ring is at r ≈ epoch/24)
    • red     →  composite (list17[i] > 0)
    • blue    →  potential twin-prime address (list17[i] == 0)
"""

import cmath
import math
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# 1.  Sieve (identical to NewGrok17.py)
# ----------------------------------------------------------------------
h = 10                                   # <-- change here for other epochs
epoch = 90 * (h * h) - 12 * h + 1        # = 8101 for h = 10
limit = epoch

# solve 90x² - 300x + (250-limit) = 0  →  positive root = new_limit
a, b, c = 90, -300, 250 - limit
disc = b**2 - 4*a*c
new_limit = (-b + cmath.sqrt(disc)) / (2*a)

list17 = [0] * (limit + 200)             # a little extra room

def drLD(x, l, m, z, lst, primitive):
    """One quadratic + arithmetic-progression marking."""
    y = 90 * (x * x) - l * x + m
    if y >= len(lst):
        return
    lst[y] += 1
    p = z + 90 * (x - 1)
    n = 1
    while True:
        pos = y + p * n
        if pos >= len(lst):
            break
        lst[pos] += 1
        n += 1

# 24 operators – exactly the calls in NewGrok17.py
for x in range(1, int(new_limit.real)):
    drLD(x,  72,  -1, 17, list17, 17)   # 17,91
    drLD(x,  72,  -1, 91, list17, 17)

    drLD(x, 108,  29, 19, list17, 17)   # 19,53
    drLD(x, 108,  29, 53, list17, 17)

    drLD(x,  72,  11, 37, list17, 17)   # 37,71
    drLD(x,  72,  11, 71, list17, 17)

    drLD(x,  18,   0, 73, list17, 17)   # 73,89
    drLD(x,  18,   0, 89, list17, 17)

    drLD(x, 102,  20, 11, list17, 17)   # 11,67
    drLD(x, 102,  20, 67, list17, 17)

    drLD(x, 138,  52, 13, list17, 17)   # 13,29
    drLD(x, 138,  52, 29, list17, 17)

    drLD(x, 102,  28, 31, list17, 17)   # 31,47
    drLD(x, 102,  28, 47, list17, 17)

    drLD(x,  48,   3, 49, list17, 17)   # 49,83
    drLD(x,  48,   3, 83, list17, 17)

    drLD(x,  78,   8, 23, list17, 17)   # 23,79
    drLD(x,  78,   8, 79, list17, 17)

    drLD(x, 132,  45,  7, list17, 17)   # 7,41
    drLD(x, 132,  45, 41, list17, 17)

    drLD(x,  78,  16, 43, list17, 17)   # 43,59
    drLD(x,  78,  16, 59, list17, 17)

    drLD(x,  42,   4, 61, list17, 17)   # 61,77
    drLD(x,  42,   4, 77, list17, 17)

# Trim the extra padding
list17 = list17[:limit]

# ----------------------------------------------------------------------
# 2.  Statistics (identical to the original program)
# ----------------------------------------------------------------------
total_marks   = sum(list17)
x_iters       = int(new_limit.real)
operators     = 24 * x_iters
unmarked      = [i for i, v in enumerate(list17) if v == 0]
quantity_holes = len(unmarked)

print(f"Epoch          : {epoch}")
print(f"Marks (sum)    : {total_marks}")
print(f"x iterations   : {x_iters}")
print(f"Operators      : {operators}")
print(f"Potential twin-prime addresses (holes) : {quantity_holes}")

# ----------------------------------------------------------------------
# 3.  Quasi-polar mapping
# ----------------------------------------------------------------------
N_RAYS = 24
ray_angle = np.linspace(0, 2*np.pi, N_RAYS, endpoint=False)   # 0, 15°, 30°, …

# radius = address / N_RAYS   →  each ring holds exactly 24 addresses
radius = np.array([i / N_RAYS for i in range(limit)], dtype=float)

# which ray does address i belong to?
ray_idx = np.arange(limit) % N_RAYS

theta = ray_angle[ray_idx]          # polar angle for every address
r     = radius                     # polar radius

# separate composites vs holes
comp_mask = np.array(list17) > 0
hole_mask = ~comp_mask

theta_comp = theta[comp_mask]
r_comp     = r[comp_mask]

theta_hole = theta[hole_mask]
r_hole     = r[hole_mask]

# ----------------------------------------------------------------------
# 4.  Plot
# ----------------------------------------------------------------------
plt.figure(figsize=(10, 10))
ax = plt.subplot(111, projection='polar')
ax.set_theta_zero_location('N')          # 0° at top
ax.set_theta_direction(-1)               # clockwise (optional, matches intuition)

# composites – small red dots
ax.scatter(theta_comp, r_comp, c='red', s=4, alpha=0.6, label='composite')

# holes – larger blue dots
ax.scatter(theta_hole, r_hole, c='dodgerblue', s=12, alpha=0.9, label='hole (potential twin)')

ax.set_rlabel_position(75)               # move radial labels a bit
ax.set_title(f"Quasi-polar sieve – h = {h} (epoch = {epoch:,})\n"
             f"{quantity_holes:,} holes, {total_marks:,} marks",
             va='bottom', pad=20)

ax.legend(loc='upper right', bbox_to_anchor=(1.15, 1.15))
plt.tight_layout()
plt.show()