#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Two Graphs:
1. Positive Space: Amplitudes (Dots + Lines + Labels)
2. Quadratic Graph: y(x) Curves + Raindrops (Lines Only) + y-dots + Red Rings at z = n
h = 20 (epoch = 35,761) — Validated with your output
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ----------------------------------------------------------------------
# 1.  Sieve – isolated lists per z
# ----------------------------------------------------------------------
h = 2                                   # Your validated case
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
    (132, 45, 7), (102, 20, 11), (138, 52, 13), 
    (72, -1, 17), (108, 29, 19), (78, 8, 23), (138, 52, 29),
    (102, 28, 31), (72, 11, 37), (132, 45, 41),
    (78, 16, 43), (102, 28, 47),
    (48, 3, 49), (108, 29, 53), (78, 16, 59),
    (42, 4, 61), (102, 20, 67), (72, 11, 71),
    (18, 0, 73), (42, 4, 77), (78, 8, 79), (48, 3, 83),
    (18, 0, 89), (72, -1, 91) 
     ]

# Store y(x) for each operator and x
y_values_per_operator = [[] for _ in range(24)]  # y_values_per_operator[k][x-1] = y for x

for k in range(24):
    l, m, z = params[k]
    z_values.append(z)
    for x in range(1, int(new_limit.real)):
        y = 90 * (x * x) - l * x + m
        if y < limit:
            y_values_per_operator[k].append(y)
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
# 3.  Positive Space Graph: Amplitudes
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

fig1 = plt.figure(figsize=(14, 11))
ax1 = fig1.add_subplot(111, projection='3d')

# Plot each ray: dots + line
colors_ray = plt.cm.tab20(np.linspace(0, 1, N_RAYS))
for k in range(N_RAYS):
    x, y, z = ray_data[k]
    ax1.scatter(x, y, z, c=[colors_ray[k]], s=8, alpha=0.7)
    ax1.plot(x, y, z, color=colors_ray[k], linewidth=1.5, alpha=0.8)

# 25th ray
ax1.scatter(x_25, y_25, z_25, c='black', s=12, alpha=0.9, label='25th: sum wave')
ax1.plot(x_25, y_25, z_25, color='black', linewidth=2, alpha=0.9)

# Reference plane
xx, yy = np.meshgrid(np.linspace(-limit/N_RAYS, limit/N_RAYS, 12),
                     np.linspace(-limit/N_RAYS, limit/N_RAYS, 12))
zz = np.zeros_like(xx)
ax1.plot_surface(xx, yy, zz, color='lightcyan', alpha=0.2)

# Red rings for holes
hole_radii = [i / N_RAYS for i in holes]
for r in hole_radii:
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    yc = r * np.sin(th)
    zc = np.zeros_like(th)
    ax1.plot(xc, yc, zc, color='red', linewidth=1.5, alpha=0.8)

# Add labels at the end of each ray
max_r = radius[-1]
for k in range(N_RAYS):
    theta_k = ray_angle[k]
    x_end = max_r * np.cos(theta_k)
    y_end = max_r * np.sin(theta_k)
    z_end = 0
    ax1.text(x_end, y_end, z_end, str(z_values[k]), color='black', fontsize=10, ha='center', va='center')

ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Amplitude')
ax1.set_title(f"Positive Space: Amplitudes – h={h}, epoch={epoch:,}\n{n_holes:,} holes → red rings")
ax1.view_init(elev=30, azim=45)
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# ----------------------------------------------------------------------
# 4.  Quadratic Graph: y(x) Curves + Raindrops (Lines Only) + y-dots + Red Rings at z = n
# ----------------------------------------------------------------------
fig2 = plt.figure(figsize=(14, 11))
ax2 = fig2.add_subplot(111, projection='3d')

max_x = min(30, x_iters)  # Show more x for curve
curve_scale = 1.0  # Scale curve height
drop_step = 1.0  # Step size for raindrops in positive y
max_drops = 8  # Max multiples per y(x)

for k in range(N_RAYS):
    theta_k = ray_angle[k]
    l, m, z = params[k]
    curve_x = []
    curve_y = []
    curve_z = []
    y_dot_x = []
    y_dot_y = []
    y_dot_z = []
    for x_idx, y_val in enumerate(y_values_per_operator[k]):
        x = x_idx + 1
        r = y_val / N_RAYS
        cx = r * np.cos(theta_k)
        cy = r * np.sin(theta_k)
        cz = curve_scale * x  # Quadratic curve in positive y
        curve_x.append(cx)
        curve_y.append(cy)
        curve_z.append(cz)
        
        # y-dot at y(x)
        y_dot_x.append(cx)
        y_dot_y.append(cy)
        y_dot_z.append(cz)
        
        # Raindrops: lines only at y + k*p in positive y direction
        p = z + 90 * (x - 1)
        pos = y_val
        drop_x = []
        drop_y = []
        drop_z = []
        for kd in range(1, max_drops + 1):
            pos += p
            if pos >= limit: break
            r_drop = pos / N_RAYS
            dx = r_drop * np.cos(theta_k)
            dy = r_drop * np.sin(theta_k)
            dz = cz + kd * drop_step  # Rise in positive y
            drop_x.append(dx)
            drop_y.append(dy)
            drop_z.append(dz)
        
        # Plot raindrops: line only
        if drop_x:
            ax2.plot(drop_x, drop_y, drop_z, color='blue', linewidth=1, alpha=0.6)
    
    # Plot quadratic curve
    if curve_x:
        ax2.plot(curve_x, curve_y, curve_z, color=colors_ray[k], linewidth=2, alpha=0.8)
    
    # Plot y-dots
    if y_dot_x:
        ax2.scatter(y_dot_x, y_dot_y, y_dot_z, c='red', s=20, alpha=0.8, marker='o')

# === Add red rings for holes at z = n (height = n) ===
for n in holes:
    r = n / N_RAYS
    z_height = n  # Height = n
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    yc = r * np.sin(th)
    zc = np.full_like(th, z_height)
    ax2.plot(xc, yc, zc, color='red', linewidth=1.5, alpha=0.8)

# Labels at curve ends
for k in range(N_RAYS):
    if y_values_per_operator[k]:
        last_y = y_values_per_operator[k][-1]
        r_end = last_y / N_RAYS
        x_end = r_end * np.cos(ray_angle[k])
        y_end = r_end * np.sin(ray_angle[k])
        z_end = curve_scale * len(y_values_per_operator[k])
        ax2.text(x_end, y_end, z_end, str(z_values[k]), color='black', fontsize=10, ha='center', va='center')

ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Positive Depth (z = n)')
ax2.set_title(f"Quadratic Graph: y(x) + Raindrops (Lines) + y-dots + Red Rings at z = n – h={h}, epoch={epoch:,}")
ax2.view_init(elev=30, azim=45)
plt.tight_layout()

plt.show()