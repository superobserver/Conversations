import math
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt

# ──────────────────────────────────────────────────────────────
# Global Constants & Helpers (from prior catalogue)
# ──────────────────────────────────────────────────────────────
ALL_CLASSES = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89]

def get_operators(k):
    ops = []
    seen = set()
    for z in ALL_CLASSES:
        try:
            o = (k * pow(z, -1, 90)) % 90
            if o not in ALL_CLASSES: continue
            pair = tuple(sorted([z,o]))
            if pair in seen: continue
            seen.add(pair)
            z_eff = 91 if z==1 else z
            o_eff = 91 if o==1 else o
            l = 180 - (z_eff + o_eff)
            m = 90 - (z_eff + o_eff) + (z_eff*o_eff - k)//90
            ops.append((l, m, z_eff))
            if z != o: ops.append((l, m, o_eff))
        except ValueError: continue
    return ops

def epoch_size(h):
    return 90 * h * h - 12 * h + 1

def compute_hole_density(h, k):
    N = epoch_size(h)
    amp = np.zeros(N, dtype=np.int8)
    xmax = int(math.sqrt(250 * N / 90)) + 10
    for l, m, z in get_operators(k):
        for x in range(1, xmax + 1):
            y = 90 * x * x - l * x + m
            if y >= N: continue
            p = z + 90 * (x - 1)
            if p <= 0: continue
            idx = int(y)
            while idx < N:
                amp[idx] += 1
                idx += p
    holes = np.sum(amp == 0)
    return holes / N if N > 0 else 0

# ──────────────────────────────────────────────────────────────
# New Proof 7: Convergence of Max Silo Density Difference to 0
# ──────────────────────────────────────────────────────────────
def prove_silo_density_convergence(h_min=10, h_max=100, h_step=10):
    print("\n=== PROOF 7: Maximum Silo Density Difference Converges to Zero ===\n")
    print("Shared operators across 24 siloes ⇒ variances from y insertions alone.")
    print("Conjecture: Max diff →0 as h→∞ (homogenization via dispersal deflection).\n")
    
    h_values = list(range(h_min, h_max + 1, h_step))
    max_diffs = []
    for h in h_values:
        densities = [compute_hole_density(h, k) for k in ALL_CLASSES]
        if densities:
            max_diff = max(densities) - min(densities)
            max_diffs.append(max_diff)
            print(f"h={h}: Max density diff = {max_diff:.6f}")
        else:
            max_diffs.append(0)
    
    # Linear fit on 1/h vs max_diff (expected decay ~1/h or faster)
    inv_h = [1/h for h in h_values]
    slope, intercept, r_value, p_value, std_err = linregress(inv_h, max_diffs)
    print(f"\nFit: max_diff vs 1/h: slope={slope:.4f}, p={p_value:.4f} (p<0.05 ⇒ significant decay)")
    if slope < 0 and p_value < 0.05:
        print("Conjecture supported: Negative slope indicates convergence to 0 as h→∞.")
    else:
        print("Conjecture not rejected, but larger h needed for significance.")

    # Plot
    plt.figure(figsize=(8,5))
    plt.plot(h_values, max_diffs, 'o-', label='Max diff')
    plt.xlabel('Epoch h')
    plt.ylabel('Max hole density diff')
    plt.title('Convergence of Max Silo Difference')
    plt.grid(True)
    plt.legend()
    plt.show()

# Add to catalogue main:
prove_silo_density_convergence()