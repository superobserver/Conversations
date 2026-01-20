import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
import multiprocessing as mp

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
    
    return np.array(list_amp)

def generate_map_grid(amplitudes, n_rows, m_cols, hole_color='#FFFF00', spiral=False):
    if spiral:
        return generate_spiral_map_grid(amplitudes, n_rows, m_cols, hole_color)
    else:
        total_pixels = n_rows * m_cols
        amps = amplitudes[:total_pixels]
        if len(amps) < total_pixels:
            amps = np.pad(amps, (0, total_pixels - len(amps)), mode='constant', constant_values=0)
        
        grid = amps.reshape(n_rows, m_cols)
        holes = (grid == 0)
        marks = grid > 0
        
        max_mark = grid[marks].max() if marks.any() else 1
        norm_marks = np.zeros_like(grid, dtype=float)
        norm_marks[marks] = grid[marks] / max_mark
        
        grid_rgb = np.zeros((n_rows, m_cols, 3), dtype=float)
        for channel in range(3):
            grid_rgb[:,:,channel][marks] = norm_marks[marks]
        
        hole_rgb = plt.cm.colors.to_rgba(hole_color)[:3]
        for channel in range(3):
            grid_rgb[:,:,channel][holes] = hole_rgb[channel]
        
        return grid_rgb

def generate_spiral_map_grid(amplitudes, n_rows, m_cols, hole_color='#FFFF00'):
    total_pixels = n_rows * m_cols
    amps = amplitudes[:total_pixels]
    if len(amps) < total_pixels:
        amps = np.pad(amps, (0, total_pixels - len(amps)), mode='constant', constant_values=0)
    
    grid = np.zeros((n_rows, m_cols))
    directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]  # right, down, left, up
    dir_idx = 0
    min_row, max_row = 0, n_rows - 1
    min_col, max_col = 0, m_cols - 1
    idx = 0
    
    while min_row <= max_row and min_col <= max_col and idx < total_pixels:
        if dir_idx == 0:  # right
            for c in range(min_col, max_col + 1):
                grid[min_row][c] = amps[idx]
                idx += 1
            min_row += 1
        elif dir_idx == 1:  # down
            for r in range(min_row, max_row + 1):
                grid[r][max_col] = amps[idx]
                idx += 1
            max_col -= 1
        elif dir_idx == 2:  # left
            for c in range(max_col, min_col - 1, -1):
                grid[max_row][c] = amps[idx]
                idx += 1
            max_row -= 1
        elif dir_idx == 3:  # up
            for r in range(max_row, min_row - 1, -1):
                grid[r][min_col] = amps[idx]
                idx += 1
            min_col += 1
        
        dir_idx = (dir_idx + 1) % 4
    
    holes = (grid == 0)
    marks = grid > 0
    
    max_mark = grid[marks].max() if marks.any() else 1
    norm_marks = np.zeros_like(grid, dtype=float)
    norm_marks[marks] = grid[marks] / max_mark
    
    grid_rgb = np.zeros((n_rows, m_cols, 3), dtype=float)
    for channel in range(3):
        grid_rgb[:,:,channel][marks] = norm_marks[marks]
    
    hole_rgb = plt.cm.colors.to_rgba(hole_color)[:3]
    for channel in range(3):
        grid_rgb[:,:,channel][holes] = hole_rgb[channel]
    
    return grid_rgb

def generate_combined_grid(amplitudes1, amplitudes2, n_rows, m_cols, spiral=False):
    if spiral:
        return generate_spiral_combined_grid(amplitudes1, amplitudes2, n_rows, m_cols)
    else:
        total_pixels = n_rows * m_cols
        amps1 = amplitudes1[:total_pixels]
        amps2 = amplitudes2[:total_pixels]
        if len(amps1) < total_pixels:
            amps1 = np.pad(amps1, (0, total_pixels - len(amps1)), mode='constant', constant_values=0)
        if len(amps2) < total_pixels:
            amps2 = np.pad(amps2, (0, total_pixels - len(amps2)), mode='constant', constant_values=0)
        
        grid1 = amps1.reshape(n_rows, m_cols)
        grid2 = amps2.reshape(n_rows, m_cols)
        
        holes1 = (grid1 == 0)
        holes2 = (grid2 == 0)
        
        grid_rgb = np.zeros((n_rows, m_cols, 3), dtype=float)
        
        color_only_11 = np.array([1.0, 0.0, 0.0])  # red
        color_only_13 = np.array([0.0, 0.0, 1.0])  # blue
        color_overlap = np.array([0.0, 1.0, 0.2])  # bright neon green
        
        only11 = holes1 & ~holes2
        only13 = ~holes1 & holes2
        both = holes1 & holes2
        
        grid_rgb[only11] = color_only_11
        grid_rgb[only13] = color_only_13
        grid_rgb[both] = color_overlap
        
        return grid_rgb

def generate_spiral_combined_grid(amplitudes1, amplitudes2, n_rows, m_cols):
    # Spiral-fill both amplitude grids separately, then compute overlaps
    grid1 = np.zeros((n_rows, m_cols))
    grid2 = np.zeros((n_rows, m_cols))
    directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]  # right, down, left, up
    dir_idx = 0
    min_row, max_row = 0, n_rows - 1
    min_col, max_col = 0, m_cols - 1
    idx = 0
    total_pixels = n_rows * m_cols
    amps1 = np.pad(amplitudes1[:total_pixels], (0, max(0, total_pixels - len(amplitudes1))), 'constant')
    amps2 = np.pad(amplitudes2[:total_pixels], (0, max(0, total_pixels - len(amplitudes2))), 'constant')
    
    while min_row <= max_row and min_col <= max_col and idx < total_pixels:
        if dir_idx == 0:  # right
            for c in range(min_col, max_col + 1):
                grid1[min_row][c] = amps1[idx]
                grid2[min_row][c] = amps2[idx]
                idx += 1
            min_row += 1
        elif dir_idx == 1:  # down
            for r in range(min_row, max_row + 1):
                grid1[r][max_col] = amps1[idx]
                grid2[r][max_col] = amps2[idx]
                idx += 1
            max_col -= 1
        elif dir_idx == 2:  # left
            for c in range(max_col, min_col - 1, -1):
                grid1[max_row][c] = amps1[idx]
                grid2[max_row][c] = amps2[idx]
                idx += 1
            max_row -= 1
        elif dir_idx == 3:  # up
            for r in range(max_row, min_row - 1, -1):
                grid1[r][min_col] = amps1[idx]
                grid2[r][min_col] = amps2[idx]
                idx += 1
            min_col += 1
        
        dir_idx = (dir_idx + 1) % 4
    
    holes1 = (grid1 == 0)
    holes2 = (grid2 == 0)
    
    grid_rgb = np.zeros((n_rows, m_cols, 3), dtype=float)
    
    color_only_11 = np.array([1.0, 0.0, 0.0])  # red
    color_only_13 = np.array([0.0, 0.0, 1.0])  # blue
    color_overlap = np.array([0.0, 1.0, 0.2])  # bright neon green
    
    only11 = holes1 & ~holes2
    only13 = ~holes1 & holes2
    both = holes1 & holes2
    
    grid_rgb[only11] = color_only_11
    grid_rgb[only13] = color_only_13
    grid_rgb[both] = color_overlap
    
    return grid_rgb

### begin generate ulam spiral ####
def generate_ulam_spiral_grid(amplitudes, n_rows, m_cols, hole_color='#FFFF00'):
    """
    Map amplitude stream into classic Ulam spiral:
    - Start at center with value amps[0] (≈ index 1 in classic Ulam)
    - Move right, then counterclockwise (up, left, down, right, ...)
    - Expanding layers
    """
    total_pixels = n_rows * m_cols
    amps = amplitudes[:total_pixels]
    if len(amps) < total_pixels:
        amps = np.pad(amps, (0, total_pixels - len(amps)), mode='constant', constant_values=0)
    
    grid = np.zeros((n_rows, m_cols))
    
    # Center starting point (for odd dimensions it's exact middle; even is biased toward top-left)
    cx = m_cols // 2
    cy = n_rows // 2
    x, y = cx, cy
    
    # Directions: right, up, left, down → counterclockwise
    directions = [(0, 1), (-1, 0), (0, -1), (1, 0)]
    dir_idx = 0
    
    # Place center
    if 0 <= x < m_cols and 0 <= y < n_rows:
        grid[y, x] = amps[0]
    
    idx = 1
    layer = 1
    while idx < total_pixels:
        # Each layer has four sides, step count increases every full layer
        for _ in range(2):  # two sides per step length
            for _ in range(layer):
                x += directions[dir_idx][1]  # column delta
                y += directions[dir_idx][0]  # row delta (note: row increases downward)
                if 0 <= x < m_cols and 0 <= y < n_rows and idx < total_pixels:
                    grid[y, x] = amps[idx]
                    idx += 1
            dir_idx = (dir_idx + 1) % 4
        
        layer += 1
    
    # Now convert to RGB visualization (same as before)
    holes = (grid == 0)
    marks = grid > 0
    
    max_mark = grid[marks].max() if marks.any() else 1
    norm_marks = np.zeros_like(grid, dtype=float)
    norm_marks[marks] = grid[marks] / max_mark
    
    grid_rgb = np.zeros((n_rows, m_cols, 3), dtype=float)
    for channel in range(3):
        grid_rgb[:,:,channel][marks] = norm_marks[marks]
    
    hole_rgb = plt.cm.colors.to_rgba(hole_color)[:3]
    for channel in range(3):
        grid_rgb[:,:,channel][holes] = hole_rgb[channel]
    
    return grid_rgb

#### end generate ulam spiral #####



def precompute_slice(args):
    classes, start, map_length = args
    return compute_amplitudes(classes, start, map_length)

def create_triple_animation(class1, class2, base_start, num_maps, map_length, n_rows, m_cols,
                            output_file='triple_animation_twins.gif', view_3d=False, spiral=False):
    # Parallel pre-computation for class1
    pool = mp.Pool(mp.cpu_count())
    args_list1 = [(class1, base_start + frame * map_length, map_length) for frame in range(num_maps)]
    all_amplitudes1 = pool.map(precompute_slice, args_list1)
    
    # Parallel pre-computation for class2
    args_list2 = [(class2, base_start + frame * map_length, map_length) for frame in range(num_maps)]
    all_amplitudes2 = pool.map(precompute_slice, args_list2)
    pool.close()
    
    if view_3d:
        fig = plt.figure(figsize=(36, 10))  # Even wider for three panels
        ax1 = fig.add_subplot(131, projection='3d')
        ax2 = fig.add_subplot(132, projection='3d')
        ax3 = fig.add_subplot(133, projection='3d')
        X, Y = np.meshgrid(range(m_cols), range(n_rows))
    else:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(36, 10))  # Triple subplots
        im1 = ax1.imshow(np.zeros((n_rows, m_cols, 3)), interpolation='nearest')
        im2 = ax2.imshow(np.zeros((n_rows, m_cols, 3)), interpolation='nearest')
        im3 = ax3.imshow(np.zeros((n_rows, m_cols, 3)), interpolation='nearest')
    
    ax1.set_title(f"Class {class1[0]} Slice 0 (Start: 0, Temp: 0.0, Holes: 0)")
    ax2.set_title(f"Class {class2[0]} Slice 0 (Start: 0, Temp: 0.0, Holes: 0)")
    ax3.set_title(f"Twins (Overlaps) Slice 0 (Start: 0, Temp: 0.0, Twins: 0)")
    ax1.set_xlabel("Column")
    ax1.set_ylabel("Row")
    ax2.set_xlabel("Column")
    ax2.set_ylabel("Row")
    ax3.set_xlabel("Column")
    ax3.set_ylabel("Row")
    ax1.axis('off')
    ax2.axis('off')
    ax3.axis('off')
    
    temperatures1 = []
    hole_counts1 = []
    temperatures2 = []
    hole_counts2 = []
    twin_temperatures = []
    twin_counts = []
    
    def animate(frame):
        amplitudes1 = all_amplitudes1[frame]
        amplitudes2 = all_amplitudes2[frame]
        
        if view_3d:
            # 3D not ideal for combined color map; fallback to 2D for simplicity
            return []
        else:
            grid_rgb1 = generate_map_grid(amplitudes1, n_rows, m_cols, spiral=spiral)
            im1.set_array(grid_rgb1)
            
            grid_rgb2 = generate_map_grid(amplitudes2, n_rows, m_cols, spiral=spiral)
            im2.set_array(grid_rgb2)
            
            grid_rgb3 = generate_combined_grid(amplitudes1, amplitudes2, n_rows, m_cols, spiral=spiral)
            im3.set_array(grid_rgb3)
        
        temp1 = np.mean(amplitudes1)
        temperatures1.append(temp1)
        hole_count1 = np.sum(amplitudes1 == 0)
        hole_counts1.append(hole_count1)
        
        temp2 = np.mean(amplitudes2)
        temperatures2.append(temp2)
        hole_count2 = np.sum(amplitudes2 == 0)
        hole_counts2.append(hole_count2)
        
        # For twins: mean amplitude where both >0, count of overlaps
        combined_amps = amplitudes1 + amplitudes2
        twin_temp = np.mean(combined_amps[ (amplitudes1 > 0) & (amplitudes2 > 0) ]) if np.any((amplitudes1 > 0) & (amplitudes2 > 0)) else 0
        twin_temperatures.append(twin_temp)
        twin_count = np.sum((amplitudes1 == 0) & (amplitudes2 == 0))
        twin_counts.append(twin_count)
        
        start_pos = base_start + frame * map_length
        ax1.set_title(f"Class {class1[0]} Slice {frame} (Start: {start_pos}, Temp: {temp1:.2f}, Holes: {hole_count1})")
        ax2.set_title(f"Class {class2[0]} Slice {frame} (Start: {start_pos}, Temp: {temp2:.2f}, Holes: {hole_count2})")
        ax3.set_title(f"Twins (Overlaps) Slice {frame} (Start: {start_pos}, Temp: {twin_temp:.2f}, Twins: {twin_count})")
        
        return [im1, im2, im3]
    
    anim = animation.FuncAnimation(fig, animate, frames=num_maps, interval=500, blit=not view_3d)
    
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Animation saved to {output_file}")
    print(f"Class {class1[0]} - Average temperatures per slice:", temperatures1)
    print(f"Class {class1[0]} - Hole counts per slice:", hole_counts1)
    print(f"Class {class2[0]} - Average temperatures per slice:", temperatures2)
    print(f"Class {class2[0]} - Hole counts per slice:", hole_counts2)
    print(f"Twin Overlaps - Average temperatures per slice:", twin_temperatures)
    print(f"Twin Overlaps - Hole counts per slice:", twin_counts)
    
    return (temperatures1, hole_counts1), (temperatures2, hole_counts2), (twin_temperatures, twin_counts)

# ================= CONFIGURABLE PARAMETERS =================
if __name__ == '__main__':
    class1 = [11]  # First class
    class2 = [13]  # Second class
    base_start = 1000000000000  # Initial starting index as requested
    num_maps = 10   # Number of successive slices
    map_length = 1000000  # Indices per map
    n_rows = 1000   # Rows in each map
    m_cols = 1000   # Columns in each map
    view_3d = False  # 3D not supported for combined color map; use False
    spiral = False  # Set to True for spiral mapping; False for standard
    create_triple_animation(class1, class2, base_start, num_maps, map_length, n_rows, m_cols, view_3d=view_3d, spiral=spiral)