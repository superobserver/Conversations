import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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

def compute_cumulative_amplitudes(classes, start_index, total_length, max_epoch):
    segment_start = start_index
    segment_end = start_index + total_length
    
    a = 90
    b = -300
    c = 250 - segment_end
    d = (b**2) - (4 * a * c)
    sol2 = (-b + math.sqrt(d)) / (2 * a) if d >= 0 else (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = int(sol2.real) + 1
    max_epoch = min(max_epoch, new_limit)  # Cap to actual limit
    
    all_ops = []
    for k in classes:
        all_ops.extend(get_operators(k))
    
    # Cumulative amplitudes per epoch frame
    cumulative_amps = np.zeros((max_epoch, total_length), dtype=int)
    
    for x in range(1, max_epoch + 1):
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
                    cumulative_amps[x-1, idx] += 1  # Add to this epoch's cumulative
                current += p
    
    # Make fully cumulative (sum prior epochs)
    for frame in range(1, max_epoch):
        cumulative_amps[frame] += cumulative_amps[frame-1]
    
    return cumulative_amps

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

def create_animation(classes, start_index, total_length, n_rows, m_cols, max_epoch, output_file='epoch_animation.gif'):
    cumulative_amps = compute_cumulative_amplitudes(classes, start_index, total_length, max_epoch)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(np.zeros((n_rows, m_cols, 3)), interpolation='nearest')
    title_text = ax.set_title("Epoch 0 (Start: 0, Temp: 0.0, Holes: 0)")
    ax.set_xlabel("Column")
    ax.set_ylabel("Row")
    ax.axis('off')
    
    temperatures = []
    hole_counts = []
    
    def animate(frame):
        amps = cumulative_amps[frame]
        grid_rgb = generate_map_grid(amps, n_rows, m_cols)
        im.set_array(grid_rgb)
        
        temp = np.mean(amps)
        temperatures.append(temp)
        
        hole_count = np.sum(amps == 0)
        hole_counts.append(hole_count)
        
        title_text.set_text(f"Epoch {frame+1} (Temp: {temp:.2f}, Holes: {hole_count})")
        return [im]
    
    anim = animation.FuncAnimation(fig, animate, frames=max_epoch, interval=200, blit=True)
    
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Animation saved to {output_file}")
    print("Average temperatures per epoch:", temperatures)
    print("Hole counts per epoch:", hole_counts)
    
    return temperatures, hole_counts

# ================= CONFIGURABLE PARAMETERS =================
classes = [11]  # Example class
start_index = 0  # Starting n
total_length = 1000000  # Fixed segment size (should be n_rows * m_cols)
n_rows = 1000
m_cols = 1000
max_epoch = 100  # Max x to iterate (controls animation length; new_limit may cap it)

# ================= RUN THE ANIMATION =================
create_animation(classes, start_index, total_length, n_rows, m_cols, max_epoch)