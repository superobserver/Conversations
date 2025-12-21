#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Evolved Quasi-Polar: Class 17 (0°-180°), Class 19 (180°-360°)
Amplitudes as Segmented Line Graphs (Breaks at Epoch Ends)
Negative z: Total amplitude for 17 (-30°), 19 (+150°), Both (+60°)
Twin Primes A224855: 90n+17 and 90n+19 both prime
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # Fixed: Import added

# ----------------------------------------------------------------------
# 1.  Sieve – isolated lists per z
# ----------------------------------------------------------------------
h = 5  # Change as needed
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
    (106, 21, 7), (110, 27, 11), (124, 40, 13), (146, 59, 17), (70, -1, 19), (74, 5, 23), (110, 33, 29), (70, 7, 31), (106, 31, 37), (110, 33, 41), (124, 40, 43), (56, 6, 47), 
    (70, 13, 49), (74, 15, 53), (110, 27, 59), (70, 13, 61), (106, 21, 67), (20, 0, 71), (34, 3, 73), (56, 6, 77),  (70, 7, 79),  (74, 5, 83), (20, 0, 89), (70, -1, 91)
    
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
# 2.  Positive Space Graph: 48 Rays (17: 0°-180°, 19: 180°-360°), Segmented by Epoch
# ----------------------------------------------------------------------
N_RAYS = 48
ray_angle = np.linspace(0, 2*np.pi, N_RAYS, endpoint=False)
radius = np.arange(limit) / N_RAYS

fig = plt.figure(figsize=(16, 12))
ax = fig.add_subplot(111, projection='3d')

# Epoch colors (cycle through for visual distinction)
epoch_colors = plt.cm.viridis(np.linspace(0, 1, int(new_limit.real)))

# Track amplitude per epoch for segmentation
amplitudes_per_epoch = [[] for _ in range(48)]  # amplitudes_per_epoch[k][x-1] = list of amplitudes for x

for k in range(48):
    theta_k = ray_angle[k]
    l, m, z = all_params[k]
    for x in range(1, int(new_limit.real)):
        y = 90 * (x * x) - l * x + m
        if y >= limit: continue
        p = z + 90 * (x - 1)
        pos = y
        amp_list = []
        while pos < limit:
            amp_list.append(all_isolated_lists[k][pos])
            pos += p
        amplitudes_per_epoch[k].append(amp_list)

# Plot segmented rays
for k in range(48):
    theta_k = ray_angle[k]
    r_start = 0
    for x_idx, amp_list in enumerate(amplitudes_per_epoch[k]):
        if not amp_list: continue
        r_vals = np.array([r_start + i for i in range(len(amp_list))]) / N_RAYS
        x_vals = r_vals * np.cos(theta_k)
        y_vals = r_vals * np.sin(theta_k)
        z_vals = np.array(amp_list)
        ax.plot(x_vals, y_vals, z_vals, color=epoch_colors[x_idx], linewidth=1.5, alpha=0.8)
        r_start += len(amp_list)

# Reference plane
xx, yy = np.meshgrid(np.linspace(-limit/N_RAYS, limit/N_RAYS, 12),
                     np.linspace(-limit/N_RAYS, limit/N_RAYS, 12))
zz = np.zeros_like(xx)
ax.plot_surface(xx, yy, zz, color='lightcyan', alpha=0.2)

# Red rings for twin holes
hole_radii = [n / N_RAYS for n in twin_holes]
for r in hole_radii:
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    yc = r * np.sin(th)
    zc = np.zeros_like(th)
    ax.plot(xc, yc, zc, color='red', linewidth=1.5, alpha=0.8)

# Labels
max_r = radius[-1]
for k in range(N_RAYS):
    theta_k = ray_angle[k]
    x_end = max_r * np.cos(theta_k)
    y_end = max_r * np.sin(theta_k)
    z_end = 0
    ax.text(x_end, y_end, z_end, str(all_z_values[k]), color='black', fontsize=9, ha='center')

# ----------------------------------------------------------------------
# 3.  Negative z: Total Amplitudes
# ----------------------------------------------------------------------
# Total amplitude for class 17
total_17 = class17_sum
r_17 = np.arange(limit) / N_RAYS
theta_17 = -np.pi / 6  # -30°
x_17 = r_17 * np.cos(theta_17)
y_17 = r_17 * np.sin(theta_17)
z_17 = -np.array(total_17)  # Negative z
ax.plot(x_17, y_17, z_17, color='blue', linewidth=2, alpha=0.9, label='Total Class 17')

# Total amplitude for class 19
total_19 = class19_sum
theta_19 = 5 * np.pi / 6  # +150° (opposite)
x_19 = r_17 * np.cos(theta_19)
y_19 = r_17 * np.sin(theta_19)
z_19 = -np.array(total_19)
ax.plot(x_19, y_19, z_19, color='red', linewidth=2, alpha=0.9, label='Total Class 19')

# Total amplitude for both
total_both = [a + b for a, b in zip(class17_sum, class19_sum)]
theta_both = np.pi / 3  # +60°
x_both = r_17 * np.cos(theta_both)
y_both = r_17 * np.sin(theta_both)
z_both = -np.array(total_both)
ax.plot(x_both, y_both, z_both, color='purple', linewidth=2, alpha=0.9, label='Total Both')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Amplitude (Negative z for Totals)')
ax.set_title(f"48 Rays: Class 17 (0°-180°), Class 19 (180°-360°) – h={h}, epoch={epoch:,}\n{n_twin_holes} twins | Epoch Segments + Negative z Totals")
ax.view_init(elev=30, azim=45)
ax.legend()
plt.tight_layout()
plt.show()