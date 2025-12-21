#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quasi-polar 3-D Sieve + Line Graphs + Ray Labels
h = 50 (epoch = 224,401) — Validated with your output
Negative Space: Quadratic y(x) curves + 'Raindrops' dots at multiples
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------------------------------------------------------------
# 1.  Sieve – isolated lists per z
# ----------------------------------------------------------------------
h = 1                                   # Your validated case
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
# 3.  3D Quasi-Polar Plot with Dots + Lines + Labels + Negative Space
# ----------------------------------------------------------------------
N_RAYS = 24
ray_angle = np.linspace(0, 2*np.pi, N_RAYS, endpoint=False)
radius = np.arange(limit) / N_RAYS

# Prepare coordinates for each ray (positive amplitude)
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

# Plot each ray: dots + line (positive z)
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

# === Negative Space: Quadratic y(x) curves + 'Raindrops' dots ===
max_x = min(20, x_iters)  # Limit x for clarity
curve_scale = 1.0  # Adjust curve height in negative z
drop_step = -1.0  # Step size for raindrops in negative z
max_drops = 5  # Max multiples per y(x)

for k in range(N_RAYS):
    theta_k = ray_angle[k]
    l, m, z = params[k]
    curve_r = []
    curve_x = []
    curve_y = []
    curve_z = []
    for x in range(1, max_x + 1):
        y = 90 * (x * x) - l * x + m
        if y >= limit: continue
        r = y / N_RAYS
        curve_r.append(r)
        cx = r * np.cos(theta_k)
        cy = r * np.sin(theta_k)
        cz = -curve_scale * x  # Quadratic curve in negative z
        curve_x.append(cx)
        curve_y.append(cy)
        curve_z.append(cz)
        
        # Raindrops: dots at y + k*p in negative z direction
        p = z + 90 * (x - 1)
        pos = y
        drop_depth = cz  # Start from curve z
        for kd in range(1, max_drops + 1):
            pos += p
            if pos >= limit: break
            r_drop = pos / N_RAYS
            dx = r_drop * np.cos(theta_k)
            dy = r_drop * np.sin(theta_k)
            dz = drop_depth - kd * drop_step  # Fall down in negative z
            ax.scatter(dx, dy, dz, c='blue', s=10, alpha=0.6)
    
    # Plot quadratic curve
    if curve_x:
        ax.plot(curve_x, curve_y, curve_z, color=colors_ray[k], linewidth=2, alpha=0.8)

# === Add labels at the end of each ray ===
max_r = radius[-1]  # End of the epoch
for k in range(N_RAYS):
    theta_k = ray_angle[k]
    x_end = max_r * np.cos(theta_k)
    y_end = max_r * np.sin(theta_k)
    z_end = 0  # Label at z=0 plane
    ax.text(x_end, y_end, z_end, str(z_values[k]), color='black', fontsize=10, ha='center', va='center')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Amplitude')
ax.set_title(f"3D Sieve – Dots + Lines + Ray Labels – h={h}, epoch={epoch:,}\n{n_holes:,} holes → red rings")
ax.view_init(elev=30, azim=45)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()