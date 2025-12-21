#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ring Geometry Graph: Class 17 (flat), Class 19 (vertical), Twin Primes (45°)
h = 5 (epoch = 2191) — Twin Primes A224855: 90n+17 and 90n+19 both prime
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
# 2.  3D Quasi-Polar Plot with 48 Rays + Squares + Red Rings at z=n
# ----------------------------------------------------------------------
N_RAYS = 48
ray_angle = np.linspace(0, 2*np.pi, N_RAYS, endpoint=False)
radius = np.arange(limit) / N_RAYS

fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111, projection='3d')

colors_ray = plt.cm.tab20(np.linspace(0, 1, N_RAYS))

# Store y(x)
y_values_per_operator = [[] for _ in range(48)]

for k in range(48):
    l, m, z = all_params[k]
    theta_k = ray_angle[k]
    curve_x = []
    curve_y = []
    curve_z = []
    for x in range(1, int(new_limit.real)):
        y_val = 90 * (x * x) - l * x + m
        if y_val >= limit: continue
        y_values_per_operator[k].append(y_val)
        
        r = y_val / N_RAYS
        cx = r * np.cos(theta_k)
        cy = r * np.sin(theta_k)
        cz = x  # z = x (epoch height)
        curve_x.append(cx)
        curve_y.append(cy)
        curve_z.append(cz)
        
        # Square
        edge = y_val
        z_height = y_val
        dx = np.cos(theta_k + np.pi/2) * edge / N_RAYS
        dy = np.sin(theta_k + np.pi/2) * edge / N_RAYS
        verts = [
            [cx, cy, z_height],
            [cx + dx, cy + dy, z_height],
            [cx + dx - np.cos(theta_k) * edge / N_RAYS, cy + dy - np.sin(theta_k) * edge / N_RAYS, z_height],
            [cx - np.cos(theta_k) * edge / N_RAYS, cy - np.sin(theta_k) * edge / N_RAYS, z_height]
        ]
        square = Poly3DCollection([verts], alpha=0.3, facecolor=colors_ray[k], edgecolor='black')
        ax.add_collection3d(square)
    
    if curve_x:
        ax.plot(curve_x, curve_y, curve_z, color=colors_ray[k], linewidth=2, alpha=0.8)

# Red rings for twin holes at z = n
for n in twin_holes:
    r = n / N_RAYS
    z_height = n
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    yc = r * np.sin(th)
    zc = np.full_like(th, z_height)
    ax.plot(xc, yc, zc, color='red', linewidth=1.5, alpha=0.8)

# Labels
for k in range(N_RAYS):
    if y_values_per_operator[k]:
        last_y = y_values_per_operator[k][-1]
        r_end = last_y / N_RAYS
        x_end = r_end * np.cos(ray_angle[k])
        y_end = r_end * np.sin(ray_angle[k])
        z_end = len(y_values_per_operator[k])
        ax.text(x_end, y_end, z_end, str(all_z_values[k]), color='black', fontsize=9, ha='center')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z = n')
ax.set_title(f"Twin Primes A224855 – h={h}, epoch={epoch:,} | {n_twin_holes} twins")
ax.view_init(elev=30, azim=45)
plt.tight_layout()

# ----------------------------------------------------------------------
# 3.  Ring Geometry Graph: Class 17 flat, Class 19 vertical, Twin Primes at 45°
# ----------------------------------------------------------------------
fig3 = plt.figure(figsize=(16, 12))
ax3 = fig3.add_subplot(111, projection='3d')

# Class 17 holes (flat, X-Y plane)
class17_holes = [n for n in range(limit) if class17_sum[n] == 0]
for n in class17_holes:
    r = n / limit * 1000  # Scale radius
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    yc = r * np.sin(th)
    zc = np.zeros_like(th)
    ax3.plot(xc, yc, zc, color='blue', linewidth=1.2, alpha=0.7)

# Class 19 holes (vertical, X-Z plane)
class19_holes = [n for n in range(limit) if class19_sum[n] == 0]
for n in class19_holes:
    r = n / limit * 1000  # Scale radius
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    zc = r * np.sin(th)
    yc = np.zeros_like(th)
    ax3.plot(xc, yc, zc, color='red', linewidth=1.2, alpha=0.7)

# Twin prime holes (at 45° relative, e.g., rotated in Y-Z plane)
for n in twin_holes:
    r = n / limit * 1000  # Scale radius
    th = np.linspace(0, 2*np.pi, 200)
    yc = r * np.cos(th)
    zc = r * np.sin(th)
    xc = np.zeros_like(th)  # In Y-Z plane at x=0
    # Rotate at 45° to make relative to others
    theta = np.pi / 4
    x_rot = xc * np.cos(theta) - yc * np.sin(theta)
    y_rot = xc * np.sin(theta) + yc * np.cos(theta)
    ax3.plot(x_rot, y_rot, zc, color='green', linewidth=1.5, alpha=0.8)

ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.set_zlabel('Z')
ax3.set_title(f"Ring Geometry: Class 17 (blue, flat), Class 19 (red, vertical), Twins (green, 45°) – h={h}, epoch={epoch:,}\nIntersect at right angles")
ax3.view_init(elev=30, azim=45)
plt.tight_layout()
plt.show()