import numpy as np
import matplotlib.pyplot as plt
from math import ceil, sqrt
import sys

# The exact 24 operators from NewGrok17.py (for 90n + 17 class)
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

def compute_epoch(h):
    return 90 * h * h - 12 * h + 1

def build_gratings(h, max_res=1024):
    epoch = compute_epoch(h)
    xmax = int((250 * epoch / 90)**0.5) + 20
    
    # Sparse storage: one set per operator
    gratings = [set() for _ in range(24)]
    combined = set()
    
    for op_idx, (l, m_val, z) in enumerate(operators):
        marked = gratings[op_idx]
        for x in range(1, xmax + 1):
            y = 90 * x * x - l * x + m_val
            if y >= epoch:
                break
            p = z + 90 * (x - 1)
            n = y
            while n < epoch:
                if n >= 0:
                    marked.add(n)
                    combined.add(n)
                n += p
                if n > epoch * 2:  # safety
                    break
    
    holes = [n for n in range(epoch) if n not in combined and n > 1000]
    
    print(f"Epoch size: {epoch:,}")
    print(f"Total marked indices: {len(combined):,}")
    print(f"Holes (primes in 90n+17): {len(holes):,} (first 10: {holes[:10]})")
    
    return gratings, combined, holes, epoch

def rasterize(marked_set, epoch, target_res=1024):
    """Turn a set of marked n into a downsampled 2D image"""
    width = target_res
    height = ceil(epoch / width)
    img = np.zeros((height, width), dtype=np.uint8)
    for n in marked_set:
        r = n // width
        c = n % width
        if r < height:
            img[r, c] = 255
    return img, width, height

def main():
    try:
        h = int(input("Enter epoch parameter h (recommended 20-300): "))
    except:
        h = 50
    
    print("Building gratings...")
    gratings, combined, holes, epoch = build_gratings(h)
    
    # Decide resolution
    if epoch < 2_000_000:  # h ≈ 150
        res = 2048
        print("High-resolution mode")
    else:
        res = 1024
        print("Downsampled mode (1024×1024)")
    
    # Combined amplitude image (downsampled)
    width = res
    height = ceil(epoch / width)
    amp_img = np.zeros((height, width), dtype=np.uint8)
    for n in combined:
        r = n // width
        c = n % width
        if r < height:
            amp_img[r, c] += 1   # will be capped later for display
    
    # Plot everything
    fig = plt.figure(figsize=(20, 18))
    fig.suptitle(f"Quadratic Sieve Diffraction Pattern — Epoch h={h} (size {epoch:,})", fontsize=16)
    
    # 24 individual gratings (6×4 grid)
    for i in range(24):
        row, col = i // 4, i % 4
        img, w, ht = rasterize(gratings[i], epoch, res)
        ax = plt.subplot(6, 4, i+1)
        ax.imshow(img, cmap='gray', origin='upper')
        ax.set_title(f"Op {i} (z={operators[i][2]})")
        ax.axis('off')
    
    # Combined amplitude
    ax_comb = plt.subplot(6, 4, 25)
    display_amp = np.clip(amp_img, 0, 10)  # cap for visibility
    ax_comb.imshow(display_amp, cmap='plasma', origin='upper')
    ax_comb.set_title("Combined (amplitude)")
    ax_comb.axis('off')
    
    # Holes as red dots on combined
    if holes:
        hole_r = [n // width for n in holes]
        hole_c = [n % width for n in holes]
        ax_comb.scatter(hole_c, hole_r, c='red', s=2, alpha=0.6)
        print(f"Plotted {len(holes):,} holes in red")
    
    plt.tight_layout()
    plt.show()
    
    # Optional: save high-res combined + holes
    plt.figure(figsize=(12,12))
    plt.imshow(display_amp, cmap='plasma', origin='upper')
    if holes:
        plt.scatter(hole_c, hole_r, c='red', s=1, alpha=0.7)
    plt.title(f"Combined Diffraction + Holes (h={h})")
    plt.axis('off')
    plt.savefig(f"sieve_diffraction_h{h}.png", dpi=300, bbox_inches='tight')
    print(f"Saved high-res image: sieve_diffraction_h{h}.png")

if __name__ == "__main__":
    main()