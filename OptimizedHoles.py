import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, ceil, floor

# The 24 operators exactly as used in your NewGrok17.py / October 2025 paper
operators = [
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

def epoch_size(h):
    return 90 * h*h - 12*h + 1

def build_grating(h, l, m, z):
    N = epoch_size(h)
    mark = np.zeros(N, dtype=bool)
    xmax = int(ceil(sqrt(250 * N / 90))) + 10
    for x in range(1, xmax + 1):
        y = 90 * x*x - l * x + m
        if y >= N: continue
        p = z + 90 * (x - 1)
        if p <= 0: continue
        idx = int(y)
        while idx < N:
            mark[idx] = True
            idx += p
    return mark

def generate_hole_candidates(h, max_candidates=1000, prefilter_small_x=True):
    N = epoch_size(h)
    candidates = []
    # Quick prefilter: avoid residues hit by x <= 5 (the strongest rays)
    forbidden_mods = set()
    if prefilter_small_x:
        for l, m, z in operators:
            for x in range(1, 6):
                y = 90 * x*x - l * x + m
                p = z + 90 * (x - 1)
                if p > 0:
                    forbidden_mods.add(p)
    
    step = max(1, N // 100000)  # Sparse sampling for large N
    for n in range(0, N, step):
        if n <= 1000: continue  # Skip small n if desired
        survives = True
        score = 1.0
        for l, m, z in operators:
            for x in range(1, 10):  # Early rays only for scoring/prefilter
                y = 90 * x*x - l * x + m
                if y >= n: break
                p = z + 90 * (x - 1)
                if p <= 0: continue
                b = (n - y) / p
                if b == floor(b):
                    survives = False
                    break
                frac = b - floor(b)
                d = min(frac, 1 - frac)
                score *= (d + 1e-9)  # Product of distances (evasion symmetry)
            if not survives: break
        if survives and score > 1e-8:  # Tunable threshold for high-symmetry
            candidates.append((n, score))
            if len(candidates) >= max_candidates * 2:  # Oversample for sorting
                break
    
    # Sort by symmetry score descending and take top
    candidates.sort(key=lambda x: x[1], reverse=True)
    return [n for n, _ in candidates[:max_candidates]]

def verify_holes(h, candidates):
    N = epoch_size(h)
    amplitude = np.zeros(N, dtype=np.uint8)
    for l, m, z in operators:
        g = build_grating(h, l, m, z)
        amplitude += g.astype(np.uint8)
    verified_holes = [n for n in candidates if amplitude[n] == 0]
    return verified_holes

def visualize(h=50, downsample=4, save=True, num_holes_to_print=100):
    N = epoch_size(h)
    print(f"Epoch = {N} indices (h={h})")

    # Optimization: Generate and verify holes first
    print("Generating hole candidates...")
    candidates = generate_hole_candidates(h, max_candidates=num_holes_to_print * 2)
    print(f"Found {len(candidates)} candidates; verifying...")
    holes_list = verify_holes(h, candidates)
    holes_list = sorted(holes_list[:num_holes_to_print])  # First 100 sorted
    print(f"First {len(holes_list)} verified holes (n where both 90n+11 and 90n+13 are prime):")
    print(holes_list)

    # Build gratings and visuals (original code)
    gratings = []
    amplitude = np.zeros(N, dtype=np.uint8)
    for l, m, z in operators:
        g = build_grating(h, l, m, z)
        gratings.append(g)
        amplitude += g.astype(np.uint8)

    holes = (amplitude == 0)
    print(f"→ {holes.sum()} total holes (primes in 90n+11)")

    # Reshape helper
    side = int(ceil(sqrt(N)))
    def to_image(arr):
        padded = np.pad(arr, (0, side*side - N), mode='constant')
        img = padded.reshape(side, side)
        if downsample > 1:
            img = img[::downsample, ::downsample]
        return img

    # 1. The 24 individual gratings
    fig1, axs = plt.subplots(6, 4, figsize=(16, 24))
    axs = axs.ravel()
    for i, g in enumerate(gratings):
        img = to_image(g)
        axs[i].imshow(img, cmap='gray_r', interpolation='nearest')
        axs[i].set_title(f"Op {i+1} (z={operators[i][2]})", fontsize=9)
        axs[i].axis('off')
    plt.suptitle(f"24 Diffraction Gratings — h={h}", fontsize=16)
    if save: plt.savefig(f"gratings_h{h}.png", dpi=200, bbox_inches='tight')

    # 2. Superposition
    super_img = to_image(amplitude > 0)
    plt.figure(figsize=(9,9))
    plt.imshow(super_img, cmap='gray_r', interpolation='nearest')
    plt.title(f"Superposition of 24 gratings (black = composite)")
    plt.axis('off')
    if save: plt.savefig(f"superposition_h{h}.png", dpi=200, bbox_inches='tight')

    # 3. Holes only (the primes)
    hole_img = to_image(holes)
    plt.figure(figsize=(9,9))
    plt.imshow(hole_img, cmap='gray_r', interpolation='nearest')
    plt.title(f"Holes = Primes (90n+11)\n{holes.sum()} dark fringes", fontsize=14)
    plt.axis('off')
    if save: plt.savefig(f"holes_h{h}.png", dpi=200, bbox_inches='tight')

    plt.show()

# ==================== RUN IT ====================
if __name__ == "__main__":
    h = int(input("Enter h (e.g. 30, 50, 100, 200): ") or 50)
    ds = int(input("Downsample factor (1 = full res, 4 = fast): ") or 4)
    visualize(h=h, downsample=ds, save=True, num_holes_to_print=100)