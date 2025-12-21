#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quasi-polar 3-D Sieve + Line Graphs (Dots + Lines)
h = 50 (epoch = 224,401) — Validated with your output
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------------------------------------------------------------
# 1.  Sieve – isolated lists per z
# ----------------------------------------------------------------------
h = 10                                   # Your validated case
epoch = 90 * (h * h) - 12 * h + 1
limit = epoch

a, b, c = 90, -300, 250 - limit
new_limit = (-b + cmath.sqrt(b**2 - 4*a*c)) / (2*a)

isolated_lists = [[0] * limit for _ in range(24)]
z_values = []

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

for k in range(24):
    l, m, z = params[k]
    z_values.append(z)
    for x in range(1, int(new_limit.real)):
        drLD(x, l, m, z, isolated_lists[k])

list_sum = [sum(isolated_lists[k][n] for k in range(24)) for n in range(limit)]

# ----------------------------------------------------------------------
# 2.  Statistics
# ----------------------------------------------------------------------
total_marks = sum(list_sum)
x_iters = int(new_limit.real)
operators = 24 * x_iters
holes = [i for i, v in enumerate(list_sum) if v == 0]
n_holes = len(holes)

print(f"Epoch: {epoch}")
print(f"Total marks: {total_marks}")
print(f"x iterations: {x_iters}")
print(f"Operators: {operators}")
print(f"Holes: {n_holes}")

# ----------------------------------------------------------------------
# 3.  3D Quasi-Polar Plot with Dots + Lines
# ----------------------------------------------------------------------
N_RAYS = 24
ray_angle = np.linspace(0, 2*np.pi, N_RAYS, endpoint=False)
radius = np.arange(limit) / N_RAYS

# Prepare coordinates for each ray
ray_data = []
for k in range(N_RAYS):
    theta_k = ray_angle[k]
    amp = np.array(isolated_lists[k])
    r = radius
    x = r * np.cos(theta_k)
    y = r * np.sin(theta_k)
    z = amp
    ray_data.append((x, y, z))

# 25th ray: aggregate sum
delta_theta = 2 * np.pi / N_RAYS / 2
theta_25 = 2 * np.pi + delta_theta
x_25 = radius * np.cos(theta_25)
y_25 = radius * np.sin(theta_25)
z_25 = np.array(list_sum)

fig = plt.figure(figsize=(14, 11))
ax = fig.add_subplot(111, projection='3d')

# Plot each ray: dots + line
colors_ray = plt.cm.tab20(np.linspace(0, 1, N_RAYS))
for k in range(N_RAYS):
    x, y, z = ray_data[k]
    # Dots
    ax.scatter(x, y, z, c=[colors_ray[k]], s=8, alpha=0.7)
    # Line (connect consecutive points)
    ax.plot(x, y, z, color=colors_ray[k], linewidth=1.5, alpha=0.8)

# 25th ray
ax.scatter(x_25, y_25, z_25, c='black', s=12, alpha=0.9, label='25th: sum wave')
ax.plot(x_25, y_25, z_25, color='black', linewidth=2, alpha=0.9)

# Reference plane
xx, yy = np.meshgrid(np.linspace(-limit/N_RAYS, limit/N_RAYS, 12),
                     np.linspace(-limit/N_RAYS, limit/N_RAYS, 12))
zz = np.zeros_like(xx)
ax.plot_surface(xx, yy, zz, color='lightcyan', alpha=0.2)

# Red rings for holes
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
ax.set_title(f"3D Sieve – Dots + Lines – h={h}, epoch={epoch:,}\n{n_holes:,} holes → red rings")
ax.view_init(elev=30, azim=45)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()