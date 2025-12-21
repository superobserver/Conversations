#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quasi-polar 3-D Sieve + Line Graphs + Ray Labels + FFT on Ray Amplitudes
h = 50 (epoch = 224,401) — Validated with your output
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fft import fft, fftfreq

# ----------------------------------------------------------------------
# 1.  Sieve – isolated lists per z
# ----------------------------------------------------------------------
h = 500                                   # Your validated case
epoch = 90 * (h * h) - 12 * h + 1
limit = epoch

a, b, c = 90, -300, 250 - limit
new_limit = (-b + cmath.sqrt(b**2 - 4*a*c)) / (2*a)

isolated_lists = [[0] * limit for _ in range(24)]
z_values = []

newvar = input("stops1")

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

for k in range(24):
    l, m, z = params[k]
    z_values.append(z)
    for x in range(1, int(new_limit.real)):
        drLD(x, l, m, z, isolated_lists[k])

list_sum = [sum(isolated_lists[k][n] for k in range(24)) for n in range(limit)]
#print(list_sum)
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
# 3.  FFT Analysis on Ray Amplitudes
# ----------------------------------------------------------------------
max_fft_n = min(limit, 200000)  # Safe for large h
n_top = 3

def compute_fft_spectrum(signal, N):
    yf = fft(signal[:N])
    xf = fftfreq(N, 1)[:N//2]
    power = 2.0 / N * np.abs(yf[0:N//2])**2
    return xf, power

spectra = []
top_harmonics = []

print("\nTop FFT harmonics (freq, power):")
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
"""
# ----------------------------------------------------------------------
# 4.  3D Quasi-Polar Plot with Dots + Lines + Labels
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
plt.show()"""
inputx = input("hmm")
# ----------------------------------------------------------------------
# 5.  Plot 2: FFT Power Spectra (24 + 1)
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

plt.suptitle(f"FFT Power Spectra – h={h}, epoch={epoch:,}", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()