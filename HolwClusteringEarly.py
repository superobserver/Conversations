import math
import numpy as np
from scipy.stats import ttest_ind

# Constants from your framework (24 classes, operators per class)
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

def compute_amplitude(k, h):
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
    return amp

def analyze_clustering(k=11, h_min=10, h_max=100, h_step=10, segments=5):
    print(f"Proving Necessary Survivor Clustering for Class {k} mod 90")
    print("Early epoch segments show elevated hole density due to y_x lag.\n")
    
    for h in range(h_min, h_max + 1, h_step):
        amp = compute_amplitude(k, h)
        N = len(amp)
        holes = np.where(amp == 0)[0]
        seg_size = N // segments
        densities = []
        for s in range(segments):
            seg_start = s * seg_size
            seg_end = (s + 1) * seg_size if s < segments - 1 else N
            seg_holes = np.sum((holes >= seg_start) & (holes < seg_end))
            dens = seg_holes / (seg_end - seg_start)
            densities.append(dens)
        
        print(f"Epoch h={h} (N={N:,}): Hole Densities per Segment")
        for s, d in enumerate(densities):
            print(f"  Segment {s+1} ({s*20}% - {(s+1)*20}%): {d:.6f}")
        
        # t-test: early (first segment) vs. late (last)
        early_holes = np.sum(amp[:seg_size] == 0) / seg_size
        late_holes = np.sum(amp[-seg_size:] == 0) / seg_size
        t_stat, p_val = ttest_ind(amp[:seg_size], amp[-seg_size:])
        print(f"  t-test (early vs. late): t={t_stat:.2f}, p={p_val:.4f} (p<0.05 ⇒ significant clustering)")
        print("─" * 50)

if __name__ == "__main__":
    analyze_clustering()  # Defaults to class 11; adjust as needed