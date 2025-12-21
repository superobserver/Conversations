#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quasi-polar 3-D visualisation – h = 10 (epoch = 8881)

* 24 rays  → each from isolated drLD markings (one list_<z> per call)
* 25th ray → total amplitude wave (sum of all list_<z>[n])
* Ring n   → address n (radius = n / 24)
* sum[n] == 0  →  draw a **red circle** through the 24 ray points at that radius
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------------------------------------------------------------
# 1.  Sieve – refactored for isolated lists per z
# ----------------------------------------------------------------------

h = input("Your EPOCH limit:")
h = int(h)

#h = 10                                   # change for other epochs
epoch = 90 * (h * h) - 12 * h + 1        # = 8881 for h = 10
limit = epoch

a, b, c = 90, -300, 250 - limit
new_limit = (-b + cmath.sqrt(b**2 - 4*a*c)) / (2*a)

# 24 isolated lists (one per drLD call, named list_<z>)
isolated_lists = [ [0] * limit for _ in range(24) ]
z_values = []  # to name them dynamically, e.g., list_17

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

# List of parameters for the 24 calls (l, m, z pairs)
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

assert len(params) == 24, "Must have exactly 24 drLD calls"

for k in range(24):
    l, m, z = params[k]
    z_values.append(z)
    for x in range(1, int(new_limit.real)):
        drLD(x, l, m, z, isolated_lists[k])

# Aggregate sum (equivalent to original list17)
list_sum = [sum(isolated_lists[k][n] for k in range(24)) for n in range(limit)]

# ----------------------------------------------------------------------
# 2.  Statistics
# ----------------------------------------------------------------------
total_marks = sum(list_sum)
x_iters     = int(new_limit.real)
operators   = 24 * x_iters
holes       = [i for i, v in enumerate(list_sum) if v == 0]
n_holes     = len(holes)

print(f"Epoch          : {epoch}")
print(f"Total marks    : {total_marks}")
print(f"x iterations   : {x_iters}")
print(f"Operators      : {operators}")
print(f"Holes (potential twins) : {n_holes}")

# Per-ray marks (for verification)
for k in range(24):
    ray_marks = sum(isolated_lists[k])
    print(f"Ray {k+1} (z={z_values[k]}) marks: {ray_marks}")

# ----------------------------------------------------------------------
# 3.  3-D quasi-polar coordinates (isolated rays)
# ----------------------------------------------------------------------
N_RAYS = 24
ray_angle = np.linspace(0, 2*np.pi, N_RAYS, endpoint=False)

# radius = address / N_RAYS
radius = np.arange(limit) / float(N_RAYS)

# For each ray k, theta_k repeated for all n, z = isolated_lists[k][n]
theta_all = np.tile(ray_angle, (limit, 1)).T  # (24, limit)
radius_all = np.tile(radius, (N_RAYS, 1))     # (24, limit)
z_all = np.array(isolated_lists)              # (24, limit)

x_all = radius_all * np.cos(theta_all)
y_all = radius_all * np.sin(theta_all)

# Flatten for scatter
x_flat = x_all.flatten()
y_flat = y_all.flatten()
z_flat = z_all.flatten()

# Colour map (blue = low/0, red = high)
z_norm = (z_flat - z_flat.min()) / (z_flat.max() - z_flat.min()) if z_flat.max() > z_flat.min() else z_flat
colors = plt.cm.coolwarm(z_norm)

# ---- 25th ray: total amplitude wave (list_sum[n] along radius) ----
delta_theta = 2 * np.pi / N_RAYS / 2.0  # half spacing offset
theta_25 = 2 * np.pi + delta_theta
x_25 = radius * np.cos(theta_25)
y_25 = radius * np.sin(theta_25)
z_25 = np.array(list_sum)  # direct sum amplitude

# ----------------------------------------------------------------------
# 4.  Plot
# ----------------------------------------------------------------------
fig = plt.figure(figsize=(13, 10))
ax  = fig.add_subplot(111, projection='3d')

# ---- scatter isolated rays -------------------------------------------------
ax.scatter(x_flat, y_flat, z_flat, c=colors, s=8, alpha=0.7)

# ---- scatter the 25th ray as a wave ------------------------------------
ax.scatter(x_25, y_25, z_25, c='black', s=12, alpha=0.9, label='25th ray: total amplitude wave')

# ---- reference plane z = 0 ---------------------------------------------
xx, yy = np.meshgrid(np.linspace(min(x_flat), max(x_flat), 12),
                     np.linspace(min(y_flat), max(y_flat), 12))
zz = np.zeros_like(xx)
ax.plot_surface(xx, yy, zz, color='lightcyan', alpha=0.2)

# ---- draw a **red ring** for **every** hole (sum[n] == 0) ------------
hole_radii = [i / N_RAYS for i in holes]
for r in hole_radii:
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    yc = r * np.sin(th)
    zc = np.zeros_like(th)
    ax.plot(xc, yc, zc, color='red', linewidth=1.5, alpha=0.8)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Amplitude')
ax.set_title(f"h = {h} – epoch = {epoch:,}\n"
             f"{n_holes:,} holes → {n_holes:,} **red rings** (one per hole)\n"
             f"24 rays: isolated drLD markings | 25th ray (black): sum wave",
             fontsize=12)

cbar = plt.colorbar(ax.scatter(x_flat, y_flat, z_flat, c=colors, s=0), ax=ax,
                    shrink=0.5, aspect=8)
cbar.set_label('Normalized amplitude')

print(f"Red rings drawn : {n_holes}")

ax.view_init(elev=30, azim=45)
ax.legend()
plt.tight_layout()
plt.show()