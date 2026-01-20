import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
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

def precompute_epoch(args):
    classes, h = args
    N_h = 90 * h**2 - 12 * h + 1
    amplitudes = compute_amplitudes(classes, 0, N_h)
    return amplitudes, N_h

def create_epoch_animation(classes, num_epochs, grid_size=100, output_file='epoch_wave_animation.gif'):
    # Precompute amplitudes cumulatively per epoch
    pool = mp.Pool(mp.cpu_count())
    args_list = [(classes, h) for h in range(1, num_epochs + 1)]
    all_data = pool.map(precompute_epoch, args_list)
    pool.close()
    
    fig, ax = plt.subplots(figsize=(10, 10))
    title = ax.set_title("Epoch 1")
    im = ax.imshow(np.zeros((grid_size, grid_size)), cmap='gray', vmin=0, vmax=1, interpolation='nearest')
    wave_lines = []
    red_scatter = ax.scatter([], [], color='red', alpha=0.5, s=10)
    
    def animate(frame):
        amplitudes, N_h = all_data[frame]
        side = math.ceil(math.sqrt(N_h))
        grid = np.zeros((side, side))
        grid.flat[:N_h] = amplitudes / (amplitudes.max() + 1e-6)  # Normalize amplitudes
        
        im.set_array(grid)
        
        # Clear previous waves and reds
        for line in wave_lines:
            line.remove()
        wave_lines.clear()
        red_scatter.set_offsets(np.empty((0,2)))
        
        # Operator transits as wave fronts
        all_ops = get_operators(classes[0])  # Assume single class for simplicity
        for l_val, m_val, z, _ in all_ops:
            for x in range(1, frame + 2):  # Up to current h
                y0 = 90 * x * x - l_val * x + m_val
                p = z + 90 * (x - 1)
                if p <= 0:
                    continue
                # Wave front as sinusoid along the transit direction
                t = np.linspace(0, 5 * p, 100)  # Several wavelengths
                wave_y = y0 + t + 0.1 * side * np.sin(2 * np.pi * t / p)  # Sinusoidal perturbation
                # Map to grid coords
                row = wave_y // side
                col = wave_y % side
                line, = ax.plot(col, row, color='blue', alpha=0.3, linewidth=0.5)
                wave_lines.append(line)
        
        # Red overlays for entanglement markings (amplitude > 0)
        marked = np.where(amplitudes > 0)[0]
        row_marked = marked // side
        col_marked = marked % side
        red_scatter.set_offsets(np.column_stack((col_marked, row_marked)))
        red_scatter.set_array(amplitudes[marked])  # Size or color by amplitude if desired
        
        ax.set_xlim(0, side)
        ax.set_ylim(0, side)
        title.set_text(f"Epoch {frame+1}, N={N_h}, Holes={np.sum(amplitudes == 0)}")
        return [im, red_scatter] + wave_lines

    anim = animation.FuncAnimation(fig, animate, frames=num_epochs, interval=500, blit=True)
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Animation saved to {output_file}")

# Configurable parameters
if __name__ == '__main__':
    classes = [11]  # Example class
    num_epochs = 200  # Number of epochs to animate (h=1 to 20)
    grid_size = 100  # Initial grid size reference (actual enlarges)
    create_epoch_animation(classes, num_epochs)