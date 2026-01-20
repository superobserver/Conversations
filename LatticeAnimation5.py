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

def generate_map_grid(amplitudes, n_rows, m_cols, hole_color='#FFFF00'):
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

def precompute_slice(args):
    classes, start, map_length = args
    return compute_amplitudes(classes, start, map_length)

def create_side_by_side_animation(class1, class2, base_start, num_maps, map_length, n_rows, m_cols,
                                 output_file='side_by_side_animation.gif', view_3d=False):
    # Parallel pre-computation for class1
    pool = mp.Pool(mp.cpu_count())
    args_list1 = [(class1, base_start + frame * map_length, map_length) for frame in range(num_maps)]
    all_amplitudes1 = pool.map(precompute_slice, args_list1)
    
    # Parallel pre-computation for class2
    args_list2 = [(class2, base_start + frame * map_length, map_length) for frame in range(num_maps)]
    all_amplitudes2 = pool.map(precompute_slice, args_list2)
    pool.close()
    
    if view_3d:
        fig = plt.figure(figsize=(24, 10))  # Wider figure for side-by-side
        ax1 = fig.add_subplot(121, projection='3d')
        ax2 = fig.add_subplot(122, projection='3d')
        X, Y = np.meshgrid(range(m_cols), range(n_rows))
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 10))  # Side-by-side subplots
        im1 = ax1.imshow(np.zeros((n_rows, m_cols, 3)), interpolation='nearest')
        im2 = ax2.imshow(np.zeros((n_rows, m_cols, 3)), interpolation='nearest')
    
    ax1.set_title(f"Class {class1[0]} Slice 0 (Start: 0, Temp: 0.0, Holes: 0)")
    ax2.set_title(f"Class {class2[0]} Slice 0 (Start: 0, Temp: 0.0, Holes: 0)")
    ax1.set_xlabel("Column")
    ax1.set_ylabel("Row")
    ax2.set_xlabel("Column")
    ax2.set_ylabel("Row")
    ax1.axis('off')
    ax2.axis('off')
    
    temperatures1 = []
    hole_counts1 = []
    temperatures2 = []
    hole_counts2 = []
    
    def animate(frame):
        amplitudes1 = all_amplitudes1[frame]
        amplitudes2 = all_amplitudes2[frame]
        
        if view_3d:
            norm_amps1 = amplitudes1 / amplitudes1.max() if amplitudes1.max() > 0 else amplitudes1
            norm_grid1 = norm_amps1.reshape(n_rows, m_cols)
            ax1.cla()
            ax1.plot_surface(X, Y, norm_grid1, cmap='viridis', vmin=0, vmax=1)
            ax1.set_zlabel('Amplitude (Height)')
            ax1.set_zlim(0, 1)
            
            norm_amps2 = amplitudes2 / amplitudes2.max() if amplitudes2.max() > 0 else amplitudes2
            norm_grid2 = norm_amps2.reshape(n_rows, m_cols)
            ax2.cla()
            ax2.plot_surface(X, Y, norm_grid2, cmap='viridis', vmin=0, vmax=1)
            ax2.set_zlabel('Amplitude (Height)')
            ax2.set_zlim(0, 1)
            return []
        else:
            grid_rgb1 = generate_map_grid(amplitudes1, n_rows, m_cols)
            im1.set_array(grid_rgb1)
            
            grid_rgb2 = generate_map_grid(amplitudes2, n_rows, m_cols)
            im2.set_array(grid_rgb2)
        
        temp1 = np.mean(amplitudes1)
        temperatures1.append(temp1)
        hole_count1 = np.sum(amplitudes1 == 0)
        hole_counts1.append(hole_count1)
        
        temp2 = np.mean(amplitudes2)
        temperatures2.append(temp2)
        hole_count2 = np.sum(amplitudes2 == 0)
        hole_counts2.append(hole_count2)
        
        start_pos = base_start + frame * map_length
        ax1.set_title(f"Class {class1[0]} Slice {frame} (Start: {start_pos}, Temp: {temp1:.2f}, Holes: {hole_count1})")
        ax2.set_title(f"Class {class2[0]} Slice {frame} (Start: {start_pos}, Temp: {temp2:.2f}, Holes: {hole_count2})")
        
        if view_3d:
            return []
        return [im1, im2]
    
    anim = animation.FuncAnimation(fig, animate, frames=num_maps, interval=500, blit=not view_3d)
    
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Animation saved to {output_file}")
    print(f"Class {class1[0]} - Average temperatures per slice:", temperatures1)
    print(f"Class {class1[0]} - Hole counts per slice:", hole_counts1)
    print(f"Class {class2[0]} - Average temperatures per slice:", temperatures2)
    print(f"Class {class2[0]} - Hole counts per slice:", hole_counts2)
    
    return (temperatures1, hole_counts1), (temperatures2, hole_counts2)

# ================= CONFIGURABLE PARAMETERS =================
if __name__ == '__main__':
    class1 = [11]  # First class
    class2 = [13]  # Second class (for side-by-side)
    base_start = 900000  # Initial starting index
    num_maps = 100   # Number of successive slices
    map_length = 10000  # Indices per map (1e6)
    n_rows = 100   # Rows in each map
    m_cols = 100   # Columns in each map
    view_3d = False  # Set True for 3D topographic view (note: 3D side-by-side may be slower)

    create_side_by_side_animation(class1, class2, base_start, num_maps, map_length, n_rows, m_cols, view_3d=view_3d)