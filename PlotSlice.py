# Updated code with slice/escape simulation
import cmath
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ALL_CLASSES = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]

def get_operators(k):
    operators = []
    seen = set()
    for z in ALL_CLASSES:
        try:
            o = (k * pow(z, -1, 90)) % 90
            if o not in ALL_CLASSES: continue
            pair = tuple(sorted([z, o]))
            if pair in seen: continue
            seen.add(pair)
            z_eff = 91 if z == 1 else z
            o_eff = 91 if o == 1 else o
            l = 180 - (z_eff + o_eff)
            m = 90 - (z_eff + o_eff) + (z_eff * o_eff - k) // 90
            operators.append((l, m, z_eff))
            if z != o:
                operators.append((l, m, o_eff))
        except ValueError:
            continue
    return operators

def compute_holes_and_slices(h, k):
    epoch = 90 * (h * h) - 12 * h + 1
    a = 90
    b = -300
    c = 250 - epoch
    d = (b**2) - (4 * a * c)
    sol2 = (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = sol2.real

    amplitude = [0] * int(epoch + 100)
    slices = []  # List of (y, pre_ratio, post_escapes)

    def drLD(x, l, m, z, listvar):
        y = 90 * (x * x) - l * x + m
        if 0 <= y < len(listvar):
            listvar[int(y)] += 1
        p = z + 90 * (x - 1)
        marked_post = 0
        for n in range(1, int((epoch - y) / p) + 1):
            yy = y + p * n
            if 0 <= yy < len(listvar):
                listvar[int(yy)] += 1
                marked_post += 1
        pre_terms = int(y)  # Finite terms < y
        post_terms = epoch - pre_terms  # "Infinite" approximation
        pre_ratio = pre_terms / epoch if epoch > 0 else 0
        post_escapes = post_terms - marked_post  # Escapes post-slice
        slices.append((y, pre_ratio, post_escapes / post_terms if post_terms > 0 else 0))

    ops = get_operators(k)
    for l, m, z in ops:
        for x in range(1, int(new_limit) + 1):
            drLD(x, l, m, z, amplitude)

    amplitude = amplitude[:int(epoch)]
    holes_indices = [i for i, amp in enumerate(amplitude) if amp == 0]
    quantity = len(holes_indices)
    density = quantity / epoch if epoch > 0 else 0
    avg_escape = np.mean([esc for _, _, esc in slices]) if slices else 0
    return h, epoch, quantity, density, slices, avg_escape

# hs = [5, 10, 20, 50]

# Narrative with slice integration
print("Chapter Addendum: Slices and Chosen Escapees")
print("Simulating for sample silo k=17, h=50: Epochs quadratic, operators slice at y, ratios b/∞, escapes as chosen primes.")
result = compute_holes_and_slices(50, 17)
print(f"Epoch: {result[1]}, Holes: {result[2]}, Density: {result[3]:.6f}")
print("Sample Slices (y, pre_ratio, post_escape_frac):", result[4][:5])  # First 5
print("Avg Escape Fraction: {result[5]:.6f} - Chosen ones evade ~this fraction per slice, compounding infinitude.")

# Plot slice distributions
plt.figure(figsize=(8, 5))
ys = [s[0] for s in result[4]]
escs = [s[2] for s in result[4]]
plt.scatter(ys, escs, c='b', label='Escape Frac per Slice')
plt.xlabel('y Slice'), plt.ylabel('Post-Slice Escape Frac'), plt.title('Chosen Escapees in Epoch')
plt.legend()
plt.savefig('slice_escapes.png')
print("Graph: 'slice_escapes.png' - Scatters show escapes post-y, finite/infinite asymmetry yielding primes.")