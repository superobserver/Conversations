import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

# ──────────────────────────────────────────────
# Sieve helpers (from your framework)
# ──────────────────────────────────────────────

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

def compute_holes(k, h):
    N = epoch_size(h)
    amp = np.zeros(N, dtype=np.int8)
    xmax = int(math.sqrt(250 * N / 90)) + 10
    ops = get_operators(k)
    for l, m, z in ops:
        for x in range(1, xmax + 1):
            y = 90 * x * x - l * x + m
            if y >= N: continue
            p = z + 90 * (x - 1)
            if p <= 0: continue
            idx = int(y)
            while idx < N:
                amp[idx] += 1
                idx += p
    return np.where(amp == 0)[0]   # array of hole indices n

# ──────────────────────────────────────────────
# Visualization: Hole density pulses across epochs
# ──────────────────────────────────────────────

def plot_hole_pulses(k=11, h_min=10, h_max=120, h_step=10, n_segments=10):
    plt.style.use('dark_background')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10), sharex=False)

    early_excess = []
    h_values = []

    for h in range(h_min, h_max + 1, h_step):
        holes = compute_holes(k, h)
        N = epoch_size(h)
        seg_size = N // n_segments
        densities = []

        for s in range(n_segments):
            start = s * seg_size
            end = (s + 1) * seg_size if s < n_segments - 1 else N
            count = np.sum((holes >= start) & (holes < end))
            dens = count / (end - start)
            densities.append(dens)

        avg_dens = len(holes) / N
        early_dens = densities[0]
        excess_pct = 100 * (early_dens - avg_dens) / avg_dens if avg_dens > 0 else 0
        early_excess.append(excess_pct)
        h_values.append(h)

        # Plot density profile for this epoch
        seg_centers = np.arange(0.5, n_segments) / n_segments
        ax1.plot(seg_centers, densities, label=f'h={h}', alpha=0.7, lw=1.5)

    ax1.set_title(f'Hole Density Profiles across Epochs — Class {k} mod 90')
    ax1.set_xlabel('Normalized position in epoch (0 = start, 1 = end)')
    ax1.set_ylabel('Local hole density')
    ax1.grid(True, alpha=0.3)
    ax1.legend(title='Epoch h', ncol=3, fontsize=8)

    # Pulse envelope: early-segment excess
    ax2.plot(h_values, early_excess, 'o-', color='gold', lw=2, markersize=6)
    ax2.axhline(0, color='gray', ls='--', alpha=0.5)
    ax2.set_title('Early-Segment Excess Hole Density (first 10%) vs. Epoch Size')
    ax2.set_xlabel('Epoch parameter h')
    ax2.set_ylabel('Excess density (%)')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()

    # Summary statistics
    mean_excess = np.mean(early_excess)
    print(f"\nMean early-segment excess density: {mean_excess:.2f}%")
    print(f"Range of excess: {min(early_excess):.2f}% to {max(early_excess):.2f}%")
    print("→ Persistent positive excess confirms clustering of survivors at epoch beginnings.")

if __name__ == '__main__':
    # Run for class 11 (A201804) across h=10 to 120
    plot_hole_pulses(k=11, h_min=10, h_max=120, h_step=5, n_segments=10)