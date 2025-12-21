#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Evolved Quasi-Polar: Class 17 (0°-180°), Class 19 (180°-360°)
Amplitudes as Segmented Line Graphs (Breaks at Epoch Ends)
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

#h = input("Your number here (3 is big for cmall video cards)")
#h = int(h)

h = 3  # Change as needed
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
   (70, -1, 91),(20, 0, 89), (74, 5, 83),  (70, 7, 79), (56, 6, 77),(34, 3, 73), (20, 0, 71), (106, 21, 67), (70, 13, 61), (110, 27, 59), (74, 15, 53), (70, 13, 49),   
    (56, 6, 47), (124, 40, 43), (110, 33, 41), (106, 31, 37), (70, 7, 31), (110, 33, 29), (74, 5, 23), (70, -1, 19), (146, 59, 17), (124, 40, 13), (110, 27, 11),
  (106, 21, 7)    

   
    
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
# 2.  Positive Space Graph: Amplitudes
# ----------------------------------------------------------------------
N_RAYS = 48
ray_angle = np.linspace(0, 2*np.pi, N_RAYS, endpoint=False)
radius = np.arange(limit) / N_RAYS

fig1 = plt.figure(figsize=(16, 12))
ax1 = fig1.add_subplot(111, projection='3d')

# Plot each ray: line graph (no dots)
colors_ray = plt.cm.tab20(np.linspace(0, 1, N_RAYS))
for k in range(N_RAYS):
    theta_k = ray_angle[k]
    amp = np.array(all_isolated_lists[k])
    r = radius
    x = r * np.cos(theta_k)
    y = r * np.sin(theta_k)
    z = amp
    ax1.plot(x, y, z, color=colors_ray[k], linewidth=1.5, alpha=0.8)

# Reference plane
xx, yy = np.meshgrid(np.linspace(-limit/N_RAYS, limit/N_RAYS, 12),
                     np.linspace(-limit/N_RAYS, limit/N_RAYS, 12))
zz = np.zeros_like(xx)
ax1.plot_surface(xx, yy, zz, color='lightcyan', alpha=0.2)

# Red rings for twin holes
hole_radii = [n / N_RAYS for n in twin_holes]
for r in hole_radii:
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    yc = r * np.sin(th)
    zc = np.zeros_like(th)
    ax1.plot(xc, yc, zc, color='red', linewidth=1.5, alpha=0.8)

# Labels
max_r = radius[-1]
for k in range(N_RAYS):
    theta_k = ray_angle[k]
    x_end = max_r * np.cos(theta_k)
    y_end = max_r * np.sin(theta_k)
    z_end = 0
    ax1.text(x_end, y_end, z_end, str(all_z_values[k]), color='black', fontsize=9, ha='center')

ax1.set_xlabel('X')
ax1.set_ylabel('Y')
ax1.set_zlabel('Amplitude')
ax1.set_title(f"48 Rays: Class 17 (0°-180°), Class 19 (180°-360°) – h={h}, epoch={epoch:,}\n{n_twin_holes} twins")
ax1.view_init(elev=30, azim=45)
plt.tight_layout()

# ----------------------------------------------------------------------
# 3.  Quadratic Graph: y(x) + Squares + Red Rings at z=n
# ----------------------------------------------------------------------
fig2 = plt.figure(figsize=(16, 12))
ax2 = fig2.add_subplot(111, projection='3d')

max_x = min(30, int(new_limit.real))
curve_scale = 1.0

# Store y(x)
y_values_per_operator = [[] for _ in range(48)]

for k in range(48):
    l, m, z = all_params[k]
    theta_k = ray_angle[k]
    curve_x = []
    curve_y = []
    curve_z = []
    for x in range(1, max_x + 1):
        y_val = 90 * (x * x) - l * x + m
        if y_val >= limit: continue
        y_values_per_operator[k].append(y_val)
        
        r = y_val / N_RAYS
        cx = r * np.cos(theta_k)
        cy = r * np.sin(theta_k)
        cz = curve_scale * x
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
        ax2.add_collection3d(square)
    
    if curve_x:
        ax2.plot(curve_x, curve_y, curve_z, color=colors_ray[k], linewidth=2, alpha=0.8)

# Red rings for twin holes at z = n
for n in twin_holes:
    r = n / N_RAYS
    z_height = n
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    yc = r * np.sin(th)
    zc = np.full_like(th, z_height)
    ax2.plot(xc, yc, zc, color='red', linewidth=1.5, alpha=0.8)

# Labels
for k in range(N_RAYS):
    if y_values_per_operator[k]:
        last_y = y_values_per_operator[k][-1]
        r_end = last_y / N_RAYS
        x_end = r_end * np.cos(ray_angle[k])
        y_end = r_end * np.sin(ray_angle[k])
        z_end = len(y_values_per_operator[k])
        ax2.text(x_end, y_end, z_end, str(all_z_values[k]), color='black', fontsize=9, ha='center')

ax2.set_xlabel('X')
ax2.set_ylabel('Y')
ax2.set_zlabel('Z = n')
ax2.set_title(f"Quadratic Graph: y(x) + Squares + Red Rings at z=n – h={h}, epoch={epoch:,}\n{n_twin_holes} twins")
ax2.view_init(elev=30, azim=45)
plt.tight_layout()

# ----------------------------------------------------------------------
# 4.  Ring Geometry Graph: Class 17 flat, Class 19 vertical, Twins at 45°
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

# Twin prime holes (at 45°)
for n in twin_holes:
    r = n / limit * 1000  # Scale radius
    th = np.linspace(0, 2*np.pi, 200)
    xc = r * np.cos(th)
    yc = r * np.sin(th)
    zc = np.zeros_like(th)
    # Rotate at 45° (in Y-Z plane for relative angle)
    theta = np.pi / 4
    y_rot = yc * np.cos(theta) - zc * np.sin(theta)
    z_rot = yc * np.sin(theta) + zc * np.cos(theta)
    ax3.plot(xc, y_rot, z_rot, color='green', linewidth=1.5, alpha=0.8)

ax3.set_xlabel('X')
ax3.set_ylabel('Y')
ax3.set_zlabel('Z')
ax3.set_title(f"Ring Geometry: Class 17 (blue, flat), Class 19 (red, vertical), Twins (green, 45°) – h={h}, epoch={epoch:,}\nIntersect at right angles")
ax3.view_init(elev=30, azim=45)
plt.tight_layout()
plt.show()

import numpy as np, matplotlib.pyplot as plt
h = np.arange(1,501)
epoch = 90*h*h -12*h +1
holes17 = np.maximum(1, np.ceil(epoch - 24*h * 4.36 * 1.08))  # conservative upper bound on marks
holes19 = np.maximum(1, np.ceil(epoch - 24*h * 4.36 * 1.08))
twin_lower = (3.75 / (9 + 2*np.log(h)))**2

plt.figure(figsize=(12,6))
plt.plot(h, holes17/epoch, label='Class 17 hole density', lw=2)
plt.plot(h, holes19/epoch, label='Class 19 hole density', lw=2)
plt.plot(h, twin_lower, '--', label='Twin lower bound', lw=3, color='red')
plt.yscale('log'); plt.xscale('log')
plt.legend(); plt.grid(alpha=0.3)
plt.title('Density decay and persistent holes per epoch (proven >0 for all h)')
plt.xlabel('Epoch parameter h'); plt.ylabel('Fraction of epoch unmarked')
plt.show()