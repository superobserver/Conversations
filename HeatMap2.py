import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def get_operators(k):
    R = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]
    operators = []
    seen_pairs = set()
    for z in R:
        try:
            o = (k * pow(z, -1, 90)) % 90
            if o not in R:
                continue
            pair = tuple(sorted([z, o]))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            z_eff = 91 if z == 1 else z
            o_eff = 91 if o == 1 else o
            l = 180 - (z_eff + o_eff)
            m = 90 - (z_eff + o_eff) + (z_eff * o_eff - k) // 90
            operators.append((l, m, z_eff, k))
            if z != o:
                operators.append((l, m, o_eff, k))
        except ValueError:
            continue
    return operators

def compute_amplitudes(classes, start_index, total_length):
    segment_start = start_index
    segment_end = start_index + total_length
    
    a = 90
    b = -300
    c = 250 - segment_end
    d = (b**2) - (4 * a * c)
    sol2 = (-b + math.sqrt(d)) / (2 * a) if d >= 0 else (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = int(sol2.real) + 1
    
    list_amp = [0] * total_length
    
    all_ops = []
    for k in classes:
        all_ops.extend(get_operators(k))
    
    for x in range(1, new_limit + 1):
        for l_val, m_val, z, primitive in all_ops:
            y0 = 90 * x * x - l_val * x + m_val
            p = z + 90 * (x - 1)
            if p <= 0 or y0 >= segment_end:
                continue
            if y0 < segment_start:
                diff = segment_start - y0
                n = (diff + p - 1) // p
                current = y0 + n * p
            else:
                current = y0
            while current < segment_end:
                if current >= segment_start:
                    idx = current - segment_start
                    list_amp[idx] += 1
                current += p
    
    return list_amp

def show_interactive_pixel_map(amplitudes, n_rows, m_cols, title="Amplitude Pixel Map"):
    total_pixels = n_rows * m_cols
    
    # Pad or truncate amplitudes to fit exactly
    amps = amplitudes[:total_pixels]
    if len(amps) < total_pixels:
        amps += [0] * (total_pixels - len(amps))  # Pad with black (holes)
    
    # Normalize to [0, 1] for colormap
    max_amp = max(amps) if max(amps) > 0 else 1
    normalized = np.array(amps) / max_amp
    
    # Reshape into n_rows x m_cols grid
    grid = normalized.reshape(n_rows, m_cols)
    
    # Create interactive figure
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(grid, cmap='viridis', interpolation='nearest', vmin=0, vmax=1)
    
    # Add colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Normalized Amplitude (0 = hole/black, 1 = brightest mark)')
    
    ax.set_title(f"{title}\n({n_rows} rows × {m_cols} columns = {total_pixels} indices)")
    ax.set_xlabel("Column index (local)")
    ax.set_ylabel("Row index (epoch-like blocks)")
    
    # Enable interactive zoom/pan
    plt.tight_layout()
    plt.show(block=True)  # Keeps window open until closed

# ================= CONFIGURABLE PARAMETERS =================
classes = [11]              # Residue class (e.g. [11] for A201804)
start_index = 0             # Starting n (index) on the line
total_length = 1000000        # Total number of consecutive indices to sieve
n_rows = 1000                # Number of rows in the grid
m_cols = 1000                # Number of columns in the grid
# total_length should equal n_rows * m_cols for perfect fit

# ================= RUN THE VISUALIZATION =================
amplitudes = compute_amplitudes(classes, start_index, total_length)
show_interactive_pixel_map(amplitudes, n_rows, m_cols,
                           title=f"Amplitude Map - Class {classes[0]} starting at n={start_index}")