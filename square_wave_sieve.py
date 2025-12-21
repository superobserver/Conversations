#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quasi-polar 3-D Sieve + Square-Wave Operator Propagation
h = 10 (epoch = 8881)
"""
import scipy.special
#print(scipy.__version__)
import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fft import fft, fftfreq
#from scipy.special import harmonic  # For H_x = 1 + 1/2 + ... + 1/x
#from scipy.special import harmonic
# ----------------------------------------------------------------------
# 1.  Sieve – isolated lists per z
# ----------------------------------------------------------------------
h = 10                                   # Change for other epochs
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
# 3.  FFT Analysis (optional, as before)
# ----------------------------------------------------------------------
# (Omit for brevity; add if needed)

# ----------------------------------------------------------------------
# 4.  3D Quasi-Polar Plot with Square Waves
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

# 25th ray
delta_theta = 2 * np.pi / N_RAYS / 2
theta_25 = 2 * np.pi + delta_theta
x_25 = radius * np.cos(theta_25)
y_25 = radius * np.sin(theta_25)
z_25 = np.array(list_sum)

fig1 = plt.figure(figsize=(13, 10))
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter(x_flat, y_flat, z_flat, c=colors, s=8, alpha=0.7)
ax1.scatter(x_25, y_25, z_25, c='black', s=12, alpha=0.9, label='25th: sum wave')

xx, yy = np.meshgrid(np.linspace(x_flat.min(), x_flat.max(), 12),
                     np.linspace(y_flat.min(), y_flat.max(), 12))
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

# Add square waves for each operator on its ray
max_wave_x = min(10, x_iters)  # Limit to avoid clutter; increase for more
colors_wave = plt.cm.viridis(np.linspace(0, 1, max_wave_x))

for k in range(N_RAYS):
    theta_k = ray_angle[k]
    for x in range(1, max_wave_x + 1):
        l, m, z = params[k]
        y = 90 * (x * x) - l * x + m
        p = z + 90 * (x - 1)
        amp = harmonic(x)  # 1 + 1/2 + ... + 1/x
        if y >= limit: continue
        
        # Square wave: steps from y onward every p
        wave_n = np.arange(y, min(y + 5 * p, limit), p)  # Short segment for vis
        wave_r = wave_n / N_RAYS
        wave_x = wave_r * np.cos(theta_k)
        wave_y = wave_r * np.sin(theta_k)
        wave_z = np.full_like(wave_r, amp)  # Constant amp for square
        
        # Step plot for square wave
        ax1.step(wave_x, wave_y, wave_z, color=colors_wave[x-1], linewidth=1.5, alpha=0.7, label=f'Ray {k+1} x={x}' if x==1 else None)

ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Amplitude')
ax1.set_title(f"3D Sieve with Square Waves – h={h}, epoch={epoch:,}\n{n_holes:,} holes → red rings")
ax1.view_init(elev=30, azim=45)
ax1.legend()
plt.tight_layout()
plt.show()