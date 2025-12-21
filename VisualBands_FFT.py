#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quasi-polar 3-D Sieve + Full FFT Power Spectra
h = 500 (epoch = 22,494,001) → validated with your output
Rays sorted by z-value: z=7 first, z=91 last
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fft import fft, fftfreq

# ----------------------------------------------------------------------
# 1.  Sieve – isolated lists per z
# ----------------------------------------------------------------------
h = 2                       # Your validated case
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

# Parameters with z-values for sorting
params_with_z = [
    (72, -1, 17), (72, -1, 91),
    (108, 29, 19), (108, 29, 53),
    (72, 11, 37), (72, 11, 71),
    (18, 0, 73), (18, 0, 89),
    (102, 20, 11), (102, 20, 67),
    (138, 52, 13), (138, 52, 29),
    (102,  28, 31), (102,  28, 47),
    (48, 3, 49), (48, 3, 83),
    (78, 8, 23), (78, 8, 79),
    (132, 45, 7), (132, 45, 41),
    (78, 16, 43), (78, 16, 59),
    (42, 4, 61), (42, 4, 77)
]

# Sort by z-value (lowest to highest)
sorted_indices = sorted(range(24), key=lambda k: params_with_z[k][2])
sorted_params = [params_with_z[i] for i in sorted_indices]

for k, (l, m, z) in enumerate(sorted_params):
    z_values.append(z)
    for x in range(1, int(new_limit.real)):
        drLD(x, l, m, z, isolated_lists[k])

# Reorder isolated_lists to match sorted z
isolated_lists = [isolated_lists[i] for i in sorted_indices]

list_sum = [sum(isolated_lists[k][n] for k in range(24)) for n in range(limit)]

# ----------------------------------------------------------------------
# 2.  Statistics
# ----------------------------------------------------------------------
total_marks = sum(list_sum)
x_iters = int(new_limit.real)
holes = [i for i, v in enumerate(list_sum) if v == 0]
n_holes = len(holes)

print(f"Epoch: {epoch}")
print(f"Total marks: {total_marks}")
print(f"x iterations: {x_iters}")
print(f"Operators: {24 * x_iters}")
print(f"Holes: {n_holes}")

# ----------------------------------------------------------------------
# 3.  FFT Analysis
# ----------------------------------------------------------------------
max_fft_n = min(limit, 200000)  # Safe for h=500
n_top = 3

def compute_fft_spectrum(signal, N):
    yf = fft(signal[:N])
    xf = fftfreq(N, 1)[:N//2]
    power = 2.0 / N * np.abs(yf[0:N//2])**2
    return xf, power

spectra = []
top_harmonics = []

print("\nTop FFT harmonics (freq, power) — sorted by z:")
for k in range(24):
    xf, power = compute_fft_spectrum(isolated_lists[k], max_fft_n)
    idx = np.argsort(power)[::-1]
    top = [(xf[i], power[i]) for i in idx[:n_top]]
    top_harmonics.append(top)
    spectra.append((xf, power))
    print(f"Ray {k+1} (z={z_values[k]}): {top}")

# Aggregate
xf_agg, power_agg = compute_fft_spectrum(list_sum, max_fft_n)
idx_agg = np.argsort(power_agg)[::-1]
top_agg = [(xf_agg[i], power_agg[i]) for i in idx_agg[:n_top]]
spectra.append((xf_agg, power_agg))
print(f"Aggregate sum wave: {top_agg}")

# ----------------------------------------------------------------------
# 4.  3D Quasi-Polar Plot
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
z_norm = (z_flat - z_flat.min()) / (z_flat.max() - z_flat.min())
colors = plt.cm.coolwarm(z_norm)

# 25th ray
delta_theta = 2 * np.pi / N_RAYS / 2
theta_25 = 2 * np.pi + delta_theta
x_25 = radius * np.cos(theta_25)
y_25 = radius * np.sin(theta_25)
z_25 = np.array(list_sum)

# ----------------------------------------------------------------------
# 5.  Plot 1: 3D Sieve
# ----------------------------------------------------------------------
fig1 = plt.figure(figsize=(13, 10))
ax1 = fig1.add_subplot(111, projection='3d')
ax1.scatter(x_flat, y_flat, z_flat, c=colors, s=6, alpha=0.7)
ax1.scatter(x_25, y_25, z_25, c='black', s=10, alpha=0.9, label='25th: sum wave')

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
    ax1.plot(xc, yc, zc, color='red', linewidth=1.2, alpha=0.8)

ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Amplitude')
ax1.set_title(f"3D Sieve – h={h}, epoch={epoch:,}\n{n_holes:,} holes → red rings")
ax1.view_init(elev=30, azim=45)
ax1.legend()
plt.tight_layout()

# ----------------------------------------------------------------------
# 6.  Plot 2: FFT Power Spectra (24 + 1) — sorted by z
# ----------------------------------------------------------------------
fig2, axes = plt.subplots(5, 5, figsize=(22, 18))
axes = axes.flatten()

for k in range(24):
    xf, power = spectra[k]
    axes[k].semilogy(xf, power, color='tab:blue', linewidth=1)
    axes[k].set_title(f"Ray {k+1} (z={z_values[k]})")
    axes[k].set_xlabel("Frequency")
    axes[k].set_ylabel("Power (log)")
    axes[k].set_xlim(0, 0.5)
    axes[k].grid(True, alpha=0.3)
    
    # Highlight 1/z, 2/z, 3/z
    z = z_values[k]
    for mult in [1, 2, 3]:
        f = mult / z
        if f <= 0.5:
            axes[k].axvline(f, color='red', linestyle='--', alpha=0.6, linewidth=1)

# Aggregate
axes[24].semilogy(xf_agg, power_agg, color='black', linewidth=1.5)
axes[24].set_title("Aggregate Sum Wave")
axes[24].set_xlabel("Frequency")
axes[24].set_ylabel("Power (log)")
axes[24].set_xlim(0, 0.5)
axes[24].grid(True, alpha=0.3)
# Highlight 1/7
axes[24].axvline(1/7, color='red', linestyle='--', alpha=0.8, linewidth=2, label='1/7')
axes[24].legend()

plt.suptitle(f"FFT Power Spectra – h={h}, epoch={epoch:,} — sorted by z", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()