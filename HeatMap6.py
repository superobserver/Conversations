import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm

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
    
    active_cohort = set()  # Unique active operators
    
    for x in range(1, new_limit + 1):
        for op in all_ops:
            l_val, m_val, z, primitive = op
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
            marked = False
            while current < segment_end:
                if current >= segment_start:
                    idx = current - segment_start
                    list_amp[idx] += 1
                    marked = True
                current += p
            if marked:
                active_cohort.add(op)  # Add unique op
    
    print(f"Intersecting operators (cohort C(segment)) for segment [{segment_start}, {segment_end}):")
    for op in sorted(active_cohort):
        print(op)
    print(f"Cohort size: {len(active_cohort)}\n")
    
    return np.array(list_amp)

def show_interactive_highlighted_map(amplitudes, n_rows, m_cols,
                                     hole_color='#FFFF00',  # Bright yellow for holes
                                     title="Highlighted Amplitude Map (Yellow = Holes)"):
    total_pixels = n_rows * m_cols
    
    # Pad or truncate
    amps = amplitudes[:total_pixels]
    if len(amps) < total_pixels:
        amps = np.pad(amps, (0, total_pixels - len(amps)), mode='constant', constant_values=0)
    
    # Create grid
    grid = amps.reshape(n_rows, m_cols)
    
    # Mask holes (0) and non-holes (>0)
    holes = (grid == 0)
    marks = grid > 0
    
    # Normalize marked amplitudes to [0,1]
    max_mark = grid[marks].max() if marks.any() else 1
    norm_marks = np.zeros_like(grid, dtype=float)
    norm_marks[marks] = grid[marks] / max_mark
    
    # Create custom colormap:
    #  - Holes (masked) → bright yellow
    #  - Marks → grayscale from black (low amp) to white (high amp)
    colors = [(0, 0, 0), (1, 1, 1)]  # black to white for marks
    mark_cmap = LinearSegmentedColormap.from_list("mark_cmap", colors)
    
    # Plot background (marks in grayscale)
    fig, ax = plt.subplots(figsize=(12, 10))
    im_marks = ax.imshow(norm_marks, cmap=mark_cmap, interpolation='nearest', vmin=0, vmax=1)
    
    # Overlay holes in bright yellow (using alpha mask)
    hole_overlay = np.zeros((n_rows, m_cols, 4))  # RGBA
    hole_overlay[holes] = plt.cm.colors.to_rgba(hole_color)  # yellow with full alpha
    ax.imshow(hole_overlay, interpolation='nearest')
    
    # Colorbar for marks (amplitude scale)
    cbar = fig.colorbar(im_marks, ax=ax, shrink=0.7)
    cbar.set_label('Normalized Amplitude (marked only)')
    
    # Title and labels
    ax.set_title(f"{title}\n({n_rows} rows × {m_cols} columns = {total_pixels} indices)\nYellow = Holes (amplitude 0)")
    ax.set_xlabel("Column index (local)")
    ax.set_ylabel("Row index (epoch-like blocks)")
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.show(block=True)  # Keeps interactive window open

# ================= CONFIGURABLE PARAMETERS =================
classes = [11]              # Example class
start_index = 900000             # Starting n
total_length = 1000000        # Total indices (should be n_rows * m_cols)
n_rows = 1000
m_cols = 1000

# ================= RUN =================
amplitudes = compute_amplitudes(classes, start_index, total_length)
show_interactive_highlighted_map(amplitudes, n_rows, m_cols,
                                 hole_color='#FFFF00',  # Change to '#00FF00' for green
                                 title=f"Class {classes[0]} - Amplitude Map with Highlighted Holes")