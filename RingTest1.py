#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quasi-polar 3-D visualisation of the quadratic twin-prime sieve
with *literal* full-zero rings drawn as red circles.

h = 10  →  epoch = 8 881
24 rays → 24 drLD operators
Ring n  →  addresses n·24 … (n+1)·24-1
A ring is drawn if **all 24 points are unmarked** (amplitude = 0).
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------------------------------------------------------------
# 1.  Sieve (identical to NewGrok17.py)
# ----------------------------------------------------------------------
h = 10                                   # change for other epochs
epoch = 90 * (h * h) - 12 * h + 1        # = 8881 for h = 10
limit = epoch

# solve 90x² - 300x + (250-limit) = 0  →  positive root = new_limit
a, b, c = 90, -300, 250 - limit
disc = b**2 - 4*a*c
new_limit = (-b + cmath.sqrt(disc)) / (2*a)

list17 = [0] * (limit + 200)             # extra safety

def drLD(x, l, m, z, lst, _):
    y = 90 * (x * x) - l * x + m
    if y >= len(lst): return
    lst[y] += 1
    p = z + 90 * (x - 1)
    n = 1
    while True:
        pos = y + p * n
        if pos >= len(lst): break
        lst[pos] += 1
        n += 1

for x in range(1, int(new_limit.real)):
    drLD(x,  72,  -1, 17, list17, 0)
    drLD(x,  72,  -1, 91, list17, 0)
    drLD(x, 108,  29, 19, list17, 0)
    drLD(x, 108,  29, 53, list17, 0)
    drLD(x,  72,  11, 37, list17, 0)
    drLD(x,  72,  11, 71, list17, 0)
    drLD(x,  18,   0, 73, list17, 0)
    drLD(x,  18,   0, 89, list17, 0)
    drLD(x, 102,  20, 11, list17, 0)
    drLD(x, 102,  20, 67, list17, 0)
    drLD(x, 138,  52, 13, list17, 0)
    drLD(x, 138,  52, 29, list17, 0)
    drLD(x, 102,  28, 31, list17, 0)
    drLD(x, 102,  28, 47, list17, 0)
    drLD(x,  48,   3, 49, list17, 0)
    drLD(x,  48,   3, 83, list17, 0)
    drLD(x,  78,   8, 23, list17, 0)
    drLD(x,  78,   8, 79, list17, 0)
    drLD(x, 132,  45,  7, list17, 0)
    drLD(x, 132,  45, 41, list17, 0)
    drLD(x,  78,  16, 43, list17, 0)
    drLD(x,  78,  16, 59, list17, 0)
    drLD(x,  42,   4, 61, list17, 0)
    drLD(x,  42,   4, 77, list17, 0)

list17 = list17[:limit]                 # trim

# ----------------------------------------------------------------------
# 2.  Basic statistics
# ----------------------------------------------------------------------
total_marks   = sum(list17)
x_iters       = int(new_limit.real)
operators     = 24 * x_iters
holes         = [i for i, v in enumerate(list17) if v == 0]
n_holes       = len(holes)

print(f"Epoch          : {epoch}")
print(f"Total marks    : {total_marks}")
print(f"x iterations   : {x_iters}")
print(f"Operators      : {operators}")
print(f"Holes (potential twins) : {n_holes}")

# ----------------------------------------------------------------------
# 3.  Ring detection (full-zero rings)
# ----------------------------------------------------------------------
N_RAYS = 24
n_rings = (limit + N_RAYS - 1) // N_RAYS          # ceiling division

zero_rings = []                                   # (ring_index, radius)
for ring in range(n_rings):
    start = ring * N_RAYS
    end   = min(start + N_RAYS, limit)
    if all(list17[i] == 0 for i in range(start, end)):
        # radius = average address of the ring / N_RAYS
        avg_addr = (start + end - 1) / 2.0
        radius   = avg_addr / N_RAYS
        zero_rings.append((ring, radius))

print(f"\nFull-zero rings found : {len(zero_rings)}")
for ring, rad in zero_rings:
    print(f"  Ring {ring:4d}  →  radius ≈ {rad:7.3f}  (addresses {ring*N_RAYS} … {min((ring+1)*N_RAYS-1, limit-1)})")

# ----------------------------------------------------------------------
# 4.  3-D quasi-polar coordinates
# ----------------------------------------------------------------------
ray_angle = np.linspace(0, 2*np.pi, N_RAYS, endpoint=False)

radius = np.arange(limit) / float(N_RAYS)          # same scaling as before
ray_idx = np.arange(limit) % N_RAYS
theta   = ray_angle[ray_idx]

x = radius * np.cos(theta)
y = radius * np.sin(theta)

# Height = marks if >0 else -1 (pits for holes)
z = np.array([v if v > 0 else -1 for v in list17])

# Colour map (blue = low/negative, red = high)
z_norm = (z - z.min()) / (z.max() - z.min()) if z.max() > z.min() else z
colors = plt.cm.coolwarm(z_norm)

# ----------------------------------------------------------------------
# 5.  Plot
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(12, 10))
ax  = fig.add_subplot(111, projection='3d')

# ---- scatter all points -------------------------------------------------
ax.scatter(x, y, z, c=colors, s=8, alpha=0.7)

# ---- reference plane z = 0 ---------------------------------------------
xx, yy = np.meshgrid(np.linspace(x.min(), x.max(), 12),
                     np.linspace(y.min(), y.max(), 12))
zz = np.zeros_like(xx)
ax.plot_surface(xx, yy, zz, color='lightcyan', alpha=0.25)

# ---- draw literal zero-rings -------------------------------------------
for ring, rad in zero_rings:
    th = np.linspace(0, 2*np.pi, 200)
    xc = rad * np.cos(th)
    yc = rad * np.sin(th)
    zc = np.zeros_like(th)
    ax.plot(xc, yc, zc, color='red', linewidth=2.5,
            label='Full-zero ring' if ring == zero_rings[0][0] else None)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Amplitude (marks)')
ax.set_title(f"h = {h}  –  epoch = {epoch:,}\n"
             f"{n_holes:,} holes  –  {len(zero_rings)} full-zero rings (red circles)",
             fontsize=12)

cbar = plt.colorbar(ax.scatter(x, y, z, c=colors, s=0), ax=ax,
                    shrink=0.5, aspect=8)
cbar.set_label('Normalized amplitude')

ax.view_init(elev=30, azim=45)
if zero_rings: ax.legend()
plt.tight_layout()
plt.show()