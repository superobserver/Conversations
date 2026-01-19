import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def get_operators(k):
    R = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89]
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
    amplitudes = {k: [0] * total_length for k in classes}
    
    for k in classes:
        segment_start = start_index
        segment_end = start_index + total_length
        
        a = 90
        b = -300
        c = 250 - segment_end
        d = (b**2) - (4 * a * c)
        sol2 = (-b + math.sqrt(d)) / (2 * a) if d >= 0 else (-b + cmath.sqrt(d)) / (2 * a)
        new_limit = int(sol2.real) + 1
        
        all_ops = get_operators(k)
        
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
                        amplitudes[k][idx] += 1
                    current += p
    
    return amplitudes

def show_interactive_overlay_map(amplitudes_dict, n_rows, m_cols, class_colors,
                                 hole_color='#FFFF00',  # Yellow for common holes/twins
                                 title="Overlay Amplitude Map (Yellow = Common Holes/Twins)"):
    total_pixels = n_rows * m_cols
    classes = list(amplitudes_dict.keys())
    
    # Prepare 2D grids for each class
    grids = {}
    for k in classes:
        amps = amplitudes_dict[k][:total_pixels]
        if len(amps) < total_pixels:
            amps += [0] * (total_pixels - len(amps))
        grid = np.array(amps).reshape(n_rows, m_cols)
        grids[k] = grid
    
    # Create RGB image (background = black)
    grid_rgb = np.zeros((n_rows, m_cols, 3), dtype=float)
    
    # Add contribution from each class
    for k in classes:
        grid = grids[k]
        marks = grid > 0
        if not marks.any():
            continue
        max_mark = grid[marks].max()
        norm_marks = np.zeros_like(grid, dtype=float)
        norm_marks[marks] = grid[marks] / max_mark
        
        color = np.array(class_colors.get(k, [1, 1, 1]))  # Default white
        for channel in range(3):
            grid_rgb[:,:,channel] += norm_marks * color[channel]
    
    # Normalize RGB so max channel doesn't exceed 1
    grid_rgb = np.clip(grid_rgb / max(1, np.max(grid_rgb)), 0, 1)
    
    # Highlight common holes (amplitude 0 in ALL classes) in yellow
    common_holes = np.all([grids[k] == 0 for k in classes], axis=0)
    hole_rgb = plt.cm.colors.to_rgba(hole_color)[:3]
    for channel in range(3):
        grid_rgb[:,:,channel][common_holes] = hole_rgb[channel]
    
    # Interactive plot
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(grid_rgb, interpolation='nearest')
    
    ax.set_title(f"{title}\n({n_rows} rows × {m_cols} columns = {total_pixels} indices)")
    ax.set_xlabel("Column index (local)")
    ax.set_ylabel("Row index (epoch-like blocks)")
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=class_colors[k], label=f"Class {k} marks") for k in classes]
    legend_elements.append(Patch(facecolor=hole_color, label="Common holes (twins)"))
    ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.show(block=True)

# ================= CONFIGURABLE PARAMETERS =================
classes = [11, 13]  # Paired classes for twins
class_colors = {11: [1, 0, 0],  # Red for class 11
                13: [0, 0, 1]}  # Blue for class 13 (overlaps = purple)
start_index = 0
total_length = 100000000
n_rows = 10000
m_cols = 10000

# ================= RUN =================
amplitudes_dict = compute_amplitudes(classes, start_index, total_length)
show_interactive_overlay_map(amplitudes_dict, n_rows, m_cols, class_colors,
                             hole_color='#FFFF00',  # Yellow for common holes
                             title=f"Twin Overlay Map - Classes {classes}")