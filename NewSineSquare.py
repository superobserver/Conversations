#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quasi-polar 3-D Sieve + Square-Wave Operator Propagation
h = 50 (epoch = 224,401) — Manual harmonic(x) — No scipy.special
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------------------------------------------------------------
# Manual harmonic function: H_x = 1 + 1/2 + ... + 1/x
# ----------------------------------------------------------------------
def harmonic(x):
    """Return the x-th harmonic number H_x = sum_{k=1}^x 1/k"""
    return sum(1.0 / k for k in range(1, x + 1))

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
# 3.  3D Quasi-Polar Plot with Square Waves
# ----------------------------------------------------------------------
N_RAYS = 24
ray_angle = np.linspace(0, 2*np.pi, N_RAYS, endpoint=False)
radius = np.arange(limit) / N_RAYS

theta_all = np.tile(ray_angle, (limit, 1)).T
radius_all = np.tile(radius, (N_RAYS, 1))
z_all = np.array(isolated_lists)

x_all = radius_all * np.cos(theta_all)
y_all = radius_all * np.sin(theta_all)

x_flat = x_all.flatten()
y_flat = y_all.flatten()
z_flat = z_all.flatten()
z_norm = (z_flat - z_flat.min()) / (z_flat.max() - z_flat.min()) if z_flat.max() > z_flat.min() else z_flat
colors = plt.cm.coolwarm(z_norm)

# 25th ray: aggregate sum
delta_theta = 2 * np.pi / N_RAYS / 2
theta_25 = 2 * np.pi + delta_theta
x_25 = radius * np.cos(theta_25)
y_25 = radius * np.sin(theta_25)
z_25 = np.array(list_sum)

fig = plt.figure(figsize=(14, 11))
ax = fig.add_subplot(111, projection='3d')

# Scatter isolated markings (dots)
ax.scatter(x_flat, y_flat, z_flat, c=colors, s=8, alpha=0.7)
ax.scatter(x_25, y_25, z_25, c='black', s=12, alpha=0.9, label='25th: sum wave')

# Reference plane
xx, yy = np.meshgrid(np.linspace(x_flat.min(), x_flat.max(), 12),
                     np.linspace(y_flat.min(), y_flat.max(), 12))
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

# === Square Waves ===
max_wave_x = min(10, x_iters)  # Show first 10 x-values to avoid clutter
colors_wave = plt.cm.viridis(np.linspace(0, 1, max_wave_x))

for k in range(N_RAYS):
    theta_k = ray_angle[k]
    l, m, z = params[k]
    for x in range(1, max_wave_x + 1):
        y = 90 * (x * x) - l * x + m
        p = z + 90 * (x - 1)
        if y >= limit: continue
        
        amp = harmonic(x)  # H_x = 1 + 1/2 + ... + 1/x
        
        # Generate wave points: start at y, step every p
        wave_n = []
        pos = y
        while pos < limit and len(wave_n) < 6:  # Limit steps for clarity
            wave_n.append(pos)
            pos += p
        if not wave_n: continue
        
        wave_r = np.array(wave_n) / N_RAYS
        wave_x = wave_r * np.cos(theta_k)
        wave_y = wave_r * np.sin(theta_k)
        wave_z = np.full_like(wave_r, amp)
        
        # Step plot for square wave
        ax.step(wave_x, wave_y, wave_z, color=colors_wave[x-1], linewidth=2, alpha=0.7,
                where='post', label=f'Ray {k+1} x={x}' if x==1 and k==0 else None)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Amplitude')
ax.set_title(f"3D Sieve + Square Waves – h={h}, epoch={epoch:,}\n"
             f"{n_holes:,} holes → red rings | Waves: H_x = harmonic(x)")
ax.view_init(elev=30, azim=45)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()