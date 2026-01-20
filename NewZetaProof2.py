import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import multiprocessing as mp

# Crystal face colors (24 distinct colors)
face_colors = plt.cm.tab20(np.linspace(0, 1, 24))

R = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]

def get_operators(k):
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
            operators.append((l, m, z_eff, k, z))  # z for face color index
            if z != o:
                operators.append((l, m, o_eff, k, o))
        except ValueError:
            continue
    return operators

def compute_amplitudes_and_traces(classes, start_index, total_length):
    segment_start = start_index
    segment_end = start_index + total_length
    
    a = 90
    b = -300
    c = 250 - segment_end
    d = (b**2) - (4 * a * c)
    sol2 = (-b + math.sqrt(d)) / (2 * a) if d >= 0 else (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = int(sol2.real) + 1
    
    list_amp = np.zeros(total_length)
    traces = []  # list of (x, y0, p, face_z) for each chain segment
    
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
            traces.append((x, y0, p, face_z))  # save for trace plotting
            while current < segment_end:
                if current >= segment_start:
                    idx = current - segment_start
                    list_amp[idx] += 1
                current += p
    
    return np.array(list_amp), traces

def precompute_epoch(args):
    classes, h = args
    N_h = 90 * h**2 - 12 * h + 1
    amplitudes, traces = compute_amplitudes_and_traces(classes, 0, N_h)
    return amplitudes, N_h, traces

def create_combined_amplitude_trace_animation(classes, num_epochs, output_file='amplitude_operator_traces.gif'):
    pool = mp.Pool(mp.cpu_count())
    args_list = [(classes, h) for h in range(1, num_epochs + 1)]
    all_data = pool.map(precompute_epoch, args_list)
    pool.close()
    
    fig, ax = plt.subplots(figsize=(12, 12))
    title = ax.set_title("Epoch 1 – Amplitude + Operator Traces")
    im = ax.imshow(np.zeros((1,1)), cmap='gray_r', vmin=0, vmax=1, interpolation='nearest', extent=[0,1,0,1])
    trace_lines = []
    red_overlay = ax.scatter([], [], color='red', alpha=0, s=12)
    
    def animate(frame):
        amplitudes, N_h, traces = all_data[frame]
        side = math.ceil(math.sqrt(N_h))
        grid = np.zeros((side, side))
        grid.flat[:N_h] = amplitudes / (amplitudes.max() + 1e-8)  # Normalize
        
        im.set_array(grid)
        im.set_extent([0, side, 0, side])
        ax.set_xlim(0, side)
        ax.set_ylim(0, side)
        
        # Clear previous
        for line in trace_lines:
            line.remove()
        trace_lines.clear()
        red_overlay.set_offsets(np.empty((0,2)))
        
        # Plot operator traces (wave fronts) colored by crystal face
        for x, y0, p, face_z in traces:
            face_idx = R.index(face_z)
            color = face_colors[face_idx]
            # Linear trace + small sinusoidal perturbation for "wave" feel
            t = np.linspace(0, min(10*p, side*1.5), 200)
            wave_y = y0 + t + 0.08 * side * np.sin(2 * np.pi * t / p)
            row = wave_y // side
            col = wave_y % side
            valid = (row >= 0) & (row < side) & (col >= 0) & (col < side)
            line, = ax.plot(col[valid], row[valid], color=color, alpha=0.35, linewidth=0.9)
            trace_lines.append(line)
        
        # Red overlay for entanglement (multi-operator positions)
        marked = np.where(amplitudes > 0)[0]
        if len(marked) > 0:
            amps_marked = amplitudes[marked]
            row_marked = marked // side
            col_marked = marked % side
            red_overlay.set_offsets(np.column_stack((col_marked, row_marked)))
            red_overlay.set_array(amps_marked / amps_marked.max())  # alpha by entanglement strength
        
        title.set_text(f"Epoch {frame+1} | N={N_h} | Holes={np.sum(amplitudes == 0)} | Operators active: {len(traces)}")
        return [im, red_overlay] + trace_lines

    anim = animation.FuncAnimation(fig, animate, frames=num_epochs, interval=800, blit=True)
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Amplitude + Operator Traces animation saved to {output_file}")

# Run parameters
if __name__ == '__main__':
    classes = [11]  # or [11,13] for twin overlay
    num_epochs = 70  # Adjust as needed (higher = more epochs, longer animation)
    create_combined_amplitude_trace_animation(classes, num_epochs)