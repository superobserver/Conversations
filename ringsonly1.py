#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ring Only Graph: Class 17 (flat) vs Class 19 (vertical), Intersect at Right Angles
h = 5 (epoch = 2191) — Twin Primes A224855: 90n+17 and 90n+19 both prime
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # Fixed: Import added

# ----------------------------------------------------------------------
# 1.  Sieve – isolated lists per z
# ----------------------------------------------------------------------
h = 1  # Change as needed
epoch = 90 * (h * h) - 12 * h + 1
limit = epoch

a, b, c = 90, -300, 250 - limit
new_limit = (-b + cmath.sqrt(b**2 - 4*a*c)) / (2*a)

# Class 17 operators (OEIS A202115: 90*k + 17 prime)
params = [
    (132, 45, 7), (102, 20, 11), (138, 52, 13), 
    (72, -1, 17), (108, 29, 19), (78, 8, 23), (138, 52, 29),
    (102, 28, 31), (72, 11, 37), (132, 45, 41),
    (78, 16, 43), (102, 28, 47),
    (48, 3, 49), (108, 29, 53), (78, 16, 59),
    (42, 4, 61), (102, 20, 67), (72, 11, 71),
    (18, 0, 73), (42, 4, 77), (78, 8, 79), (48, 3, 83),
    (18, 0, 89), (72, -1, 91) 
]

# Class 19 operators (OEIS A196000: 90*k + 19 prime) — Your corrected list
params_19 = [
    (70, -1, 19), (70, -1, 91),
    (106, 31, 37), (34, 3, 73),
    (74, 15, 53), (146, 59, 17),
    (110, 27, 11), (110, 27, 59),
    (110, 33, 29), (110, 33, 41),
    (56, 6, 47), (56, 6, 77),
    (74, 5, 23), (74, 5, 83),
    (124, 40, 13), (124, 40, 43),
    (70, 7, 31), (70, 7, 79),
    (70, 13, 49), (70, 13, 61),
    (106, 21, 7), (106, 21, 67),
    (20, 0, 71), (20, 0, 89)
]

# Combine into 48 operators
all_params = params + params_19
all_z_values = [z for l, m, z in all_params]

# Isolated amplitude lists
all_isolated_lists = [[0] * limit for _ in range(48)]

def drLD(x, l, m, z, lst):
    y = 90 * (x * x) - l * x + m
    if y >= len(lst): return
    lst[int(y)] += 1
    p = z + 90 * (x - 1)
    n = 1
    while True:
        pos = y + p * n
        if pos >= len(lst): break
        lst[int(pos)] += 1
        n += 1

# Run sieve
for k in range(48):
    l, m, z = all_params[k]
    for x in range(1, int(new_limit.real)):
        drLD(x, l, m, z, all_isolated_lists[k])

# Class sums
class17_sum = [sum(all_isolated_lists[k][n] for k in range(24)) for n in range(limit)]
class19_sum = [sum(all_isolated_lists[k][n] for k in range(24, 48)) for n in range(limit)]

# Twin prime holes (common unmarked)
twin_holes = [n for n in range(limit) if class17_sum[n] == 0 and class19_sum[n] == 0]
n_twin_holes = len(twin_holes)

print(f"Epoch: {epoch}")
print(f"Class 17 holes: {sum(1 for v in class17_sum if v == 0)}")
print(f"Class 19 holes: {sum(1 for v in class19_sum if v == 0)}")
print(f"Twin primes (A224855): {n_twin_holes}")

# ----------------------------------------------------------------------
# 2.  Ring Geometry Graph: Class 17 flat, Class 19 vertical
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111, projection='3d')

# Class 17 holes (flat, X-Y plane)
class17_holes = [n for n in range(limit) if class17_sum[n] == 0]
for n in class17_holes:
    r = n / limit * 1000  # Scale radius
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    yc = r * np.sin(th)
    zc = np.zeros_like(th)
    ax.plot(xc, yc, zc, color='blue', linewidth=1.2, alpha=0.7)

# Class 19 holes (vertical, X-Z plane)
class19_holes = [n for n in range(limit) if class19_sum[n] == 0]
for n in class19_holes:
    r = n / limit * 1000  # Scale radius
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    zc = r * np.sin(th)
    yc = np.zeros_like(th)
    ax.plot(xc, yc, zc, color='red', linewidth=1.2, alpha=0.7)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title(f"Ring Geometry: Class 17 (flat X-Y) vs Class 19 (vertical X-Z) – h={h}, epoch={epoch:,}\nIntersect at 90°")
ax.view_init(elev=30, azim=45)
plt.tight_layout()
plt.show()