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

def create_animation(classes, base_start, num_maps, map_length, n_rows, m_cols,
                     output_file='amplitude_time_series.gif', view_3d=False):
    # Parallel pre-computation
    pool = mp.Pool(mp.cpu_count())
    args_list = [(classes, base_start + frame * map_length, map_length) for frame in range(num_maps)]
    all_amplitudes = pool.map(precompute_slice, args_list)
    pool.close()
    
    if view_3d:
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        X, Y = np.meshgrid(range(m_cols), range(n_rows))
    else:
        fig, ax = plt.subplots(figsize=(12, 10))
        im = ax.imshow(np.zeros((n_rows, m_cols, 3)), interpolation='nearest')
    
    title_text = ax.set_title("Slice 0 (Start: 0, Temp: 0.0, Holes: 0)")
    ax.set_xlabel("Column")
    ax.set_ylabel("Row")
    ax.axis('off')
    
    temperatures = []
    hole_counts = []
    
    def animate(frame):
        amplitudes = all_amplitudes[frame]
        if view_3d:
            norm_amps = amplitudes / amplitudes.max() if amplitudes.max() > 0 else amplitudes
            norm_grid = norm_amps.reshape(n_rows, m_cols)
            ax.cla()
            ax.plot_surface(X, Y, norm_grid, cmap='viridis', vmin=0, vmax=1)
            ax.set_zlabel('Amplitude (Height)')
            ax.set_zlim(0, 1)
        else:
            grid_rgb = generate_map_grid(amplitudes, n_rows, m_cols)
            im.set_array(grid_rgb)
        
        temp = np.mean(amplitudes)
        temperatures.append(temp)
        
        hole_count = np.sum(amplitudes == 0)
        hole_counts.append(hole_count)
        
        title_text.set_text(f"Slice {frame} (Start: {base_start + frame * map_length}, Temp: {temp:.2f}, Holes: {hole_count})")
        if view_3d:
            return []
        return [im]
    
    anim = animation.FuncAnimation(fig, animate, frames=num_maps, interval=500, blit=not view_3d)
    
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Animation saved to {output_file}")
    print("Average temperatures per slice:", temperatures)
    print("Hole counts per slice:", hole_counts)
    
    return temperatures, hole_counts

# ================= CONFIGURABLE PARAMETERS =================
if __name__ == '__main__':
    classes = [11]  # Example class
    base_start = 900000000000  # Initial starting index
    num_maps = 300   # Number of successive slices
    map_length = 1000000  # Indices per map (1e6)
    n_rows = 1000   # Rows in each map
    m_cols = 1000   # Columns in each map
    view_3d = False  # Set True for 3D topographic view

    create_animation(classes, base_start, num_maps, map_length, n_rows, m_cols, view_3d=view_3d)
