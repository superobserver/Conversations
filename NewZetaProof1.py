import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import multiprocessing as mp

# Crystal face colors (24 distinct colors for the 24 residues)
face_colors = plt.cm.tab20(np.linspace(0, 1, 24))

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
            operators.append((l, m, z_eff, k, z))  # add z for face color
            if z != o:
                operators.append((l, m, o_eff, k, o))
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
        for l_val, m_val, z, primitive, face_z in all_ops:
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
    
    return np.array(list_amp), all_ops

def precompute_epoch(args):
    classes, h = args
    N_h = 90 * h**2 - 12 * h + 1
    amplitudes, ops = compute_amplitudes(classes, 0, N_h)
    return amplitudes, N_h, ops

def create_epoch_crystal_animation(classes, num_epochs, output_file='epoch_crystal_illumination.gif'):
    R = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]
    pool = mp.Pool(mp.cpu_count())
    args_list = [(classes, h) for h in range(1, num_epochs + 1)]
    all_data = pool.map(precompute_epoch, args_list)
    pool.close()
    
    fig, ax = plt.subplots(figsize=(12, 12))
    title = ax.set_title("Epoch 1 – Crystal Illumination")
    im = ax.imshow(np.zeros((1,1)), cmap='gray', vmin=0, vmax=1, interpolation='nearest')
    wave_lines = []
    red_overlay = ax.scatter([], [], color='red', alpha=0, s=10)
    
    # Crystal face legend (inset)
    legend_ax = fig.add_axes([0.02, 0.02, 0.2, 0.2])
    legend_ax.axis('off')
    for i, z in enumerate(R):
        legend_ax.add_patch(plt.Rectangle((0, i/24), 0.2, 1/24, color=face_colors[i]))
        legend_ax.text(0.25, i/24 + 0.02/24, str(z), va='center', fontsize=8)
    legend_ax.set_title('Zeta-Crystal Faces', fontsize=9)
    
    def animate(frame):
        amplitudes, N_h, ops = all_data[frame]
        side = math.ceil(math.sqrt(N_h))
        grid = np.zeros((side, side))
        grid.flat[:N_h] = amplitudes / (amplitudes.max() + 1e-6)
        
        im.set_array(grid)
        ax.set_xlim(0, side)
        ax.set_ylim(0, side)
        
        # Clear previous
        for line in wave_lines:
            line.remove()
        wave_lines.clear()
        red_overlay.set_offsets(np.empty((0,2)))
        
        # Wave fronts colored by crystal face
        for l_val, m_val, z, primitive, face_z in ops:
            face_idx = R.index(face_z)
            color = face_colors[face_idx]
            for x in range(1, frame + 2):
                y0 = 90 * x * x - l_val * x + m_val
                p = z + 90 * (x - 1)
                if p <= 0:
                    continue
                t = np.linspace(0, 6 * p, 150)
                wave_y = y0 + t + 0.08 * side * np.sin(2 * np.pi * t / p)
                row = wave_y // side
                col = wave_y % side
                line, = ax.plot(col, row, color=color, alpha=0.4, linewidth=0.8)
                wave_lines.append(line)
        
        # Red overlay for multi-face entanglement
        marked = np.where(amplitudes > 0)[0]
        if len(marked) > 0:
            amps_marked = amplitudes[marked]
            row_marked = marked // side
            col_marked = marked % side
            red_overlay.set_offsets(np.column_stack((col_marked, row_marked)))
            red_overlay.set_array(amps_marked / amps_marked.max())  # alpha by amplitude count
        
        title.set_text(f"Epoch {frame+1} | N={N_h} | Holes={np.sum(amplitudes == 0)}")
        return [im, red_overlay] + wave_lines

    anim = animation.FuncAnimation(fig, animate, frames=num_epochs, interval=600, blit=True)
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Crystal-illuminated epoch animation saved to {output_file}")

# Run parameters
if __name__ == '__main__':
    classes = [11]  # Can add [11,13] for twin view
    num_epochs = 30  # Adjust for longer animation
    create_epoch_crystal_animation(classes, num_epochs)