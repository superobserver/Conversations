import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Constants from your framework (24 operators for class 11, as example; extend as needed)
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
    return 90 * h * h - 12 * h + 1

def build_stitch_pattern(h, l, m, z):
    N = epoch_size(h)
    stitch = np.zeros(N, dtype=bool)
    xmax = int(math.sqrt(250 * N / 90)) + 10
    for x in range(1, xmax + 1):
        y = 90 * x * x - l * x + m
        if y >= N: continue
        p = z + 90 * (x - 1)
        if p <= 0: continue
        idx = int(y)
        while idx < N:
            stitch[idx] = True
            idx += p
    return stitch

def stitch_embroidery(h=50, save=True):
    N = epoch_size(h)
    print(f"Epoch = {N} indices (h={h})")

    # Build 24 individual stitch patterns + summed amplitude
    patterns = []
    amplitude = np.zeros(N, dtype=np.uint8)
    for i, (l, m, z) in enumerate(operators):
        pat = build_stitch_pattern(h, l, m, z)
        patterns.append(pat)
        amplitude += pat.astype(np.uint8)

    prime_holes = np.where(amplitude == 0)[0]
    print(f"→ {len(prime_holes)} prime holes (diamonds)")

    # Reshape to n-by-n grid (side ≈ sqrt(N))
    side = int(math.ceil(math.sqrt(N)))
    def to_grid(arr):
        padded = np.pad(arr, (0, side*side - N), mode='constant')
        return padded.reshape(side, side)

    # 1. Individual operator stitches (24 patterns)
    fig1, axs = plt.subplots(6, 4, figsize=(16, 24))
    axs = axs.ravel()
    for i, pat in enumerate(patterns):
        grid = to_grid(pat)
        axs[i].imshow(grid, cmap='binary', interpolation='nearest')
        axs[i].set_title(f"Op {i+1} (z={operators[i][2]})", fontsize=9)
        axs[i].axis('off')
    plt.suptitle(f"24 Operator Stitch Patterns — h={h}", fontsize=16)
    if save: plt.savefig(f"operator_stitches_h{h}.png", dpi=200, bbox_inches='tight')

    # 2. Sum-stitch (overlaid patterns, opacity ~ amplitude)
    sum_grid = to_grid(amplitude / np.max(amplitude))
    plt.figure(figsize=(9,9))
    plt.imshow(sum_grid, cmap='Greys', interpolation='nearest')
    plt.title(f"Sum-Stitch of 24 Operators (darker = higher amplitude)")
    plt.axis('off')
    if save: plt.savefig(f"sum_stitch_h{h}.png", dpi=200, bbox_inches='tight')

    # 3. Prime holes with diamonds
    hole_grid = to_grid(amplitude == 0)
    plt.figure(figsize=(9,9))
    plt.imshow(hole_grid, cmap='gray_r', interpolation='nearest')
    for hole in prime_holes[:100]:  # Limit to first 100 for clarity
        row, col = divmod(hole, side)
        diamond = Polygon([[col-0.25, row], [col, row-0.25], [col+0.25, row], [col, row+0.25]], closed=True, fill=True, color='red')
        plt.gca().add_patch(diamond)
    plt.title(f"Prime Holes Marked with Diamonds (red) — {len(prime_holes)} voids")
    plt.axis('off')
    if save: plt.savefig(f"prime_diamonds_h{h}.png", dpi=200, bbox_inches='tight')

    plt.show()

if __name__ == "__main__":
    h = int(input("Enter h (e.g. 30, 50): ") or 50)
    stitch_embroidery(h=h, save=True)