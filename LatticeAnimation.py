import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
    
    return np.array(list_amp)

def generate_map_grid(amplitudes, n_rows, m_cols, hole_color='#FFFF00'):
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
    
    # Create RGB image
    grid_rgb = np.zeros((n_rows, m_cols, 3), dtype=float)
    
    # Grayscale for marks (0 black to 1 white)
    for channel in range(3):
        grid_rgb[:,:,channel][marks] = norm_marks[marks]
    
    # Highlight holes in yellow
    hole_rgb = plt.cm.colors.to_rgba(hole_color)[:3]
    for channel in range(3):
        grid_rgb[:,:,channel][holes] = hole_rgb[channel]
    
    return grid_rgb

def create_animation(classes, base_start, num_maps, map_length, n_rows, m_cols, output_file='amplitude_time_series.gif'):
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(np.zeros((n_rows, m_cols, 3)), interpolation='nearest')
    title_text = ax.set_title("Amplitude Map - Slice 0 (Start: 0, Temp: 0.0)")
    ax.set_xlabel("Column")
    ax.set_ylabel("Row")
    ax.axis('off')
    
    temperatures = []  # To store average amplitudes
    
    def animate(frame):
        start = base_start + frame * map_length
        amplitudes = compute_amplitudes(classes, start, map_length)
        grid_rgb = generate_map_grid(amplitudes, n_rows, m_cols)
        
        temp = np.mean(amplitudes)
        temperatures.append(temp)
        
        im.set_array(grid_rgb)
        title_text.set_text(f"Amplitude Map - Slice {frame} (Start: {start}, Temp: {temp:.2f})")
        return [im]
    
    anim = animation.FuncAnimation(fig, animate, frames=num_maps, interval=500, blit=True)
    
    # Save as GIF
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Animation saved to {output_file}")
    print("Average temperatures per slice:", temperatures)
    
    return temperatures

# ================= CONFIGURABLE PARAMETERS =================
classes = [11]  # Example class
base_start = 0  # Initial starting index
num_maps = 20   # Number of successive maps (slices) to generate
map_length = 1000000  # Indices per map (1e6)
n_rows = 1000   # Rows in each map
m_cols = 1000   # Columns in each map (map_length must = n_rows * m_cols)

# ================= RUN THE ANIMATION =================
create_animation(classes, base_start, num_maps, map_length, n_rows, m_cols)