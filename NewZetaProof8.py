import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import multiprocessing as mp
from matplotlib.patches import Rectangle
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
            operators.append((l, m, z_eff, k, z))  # z for face color
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
    traces = []  # (x, y0, p, face_z)
    
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
            traces.append((x, y0, p, face_z))
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


# ... (keep all previous imports and functions: get_operators, compute_amplitudes_and_traces, precompute_epoch)

def create_dual_view_animation(classes, num_epochs, show_traces_panel=True, show_red_intensity=True, output_file='dual_view_epoch_animation_centered.gif'):
    pool = mp.Pool(mp.cpu_count())
    args_list = [(classes, h) for h in range(1, num_epochs + 1)]
    all_data = pool.map(precompute_epoch, args_list)
    pool.close()
    
    if show_traces_panel:
        fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(18, 9))
    else:
        fig, ax_left = plt.subplots(figsize=(12, 12))
        ax_right = None
    
    title_left = ax_left.set_title("Epoch 1 – Amplitude Field")
    im = ax_left.imshow(np.zeros((1,1)), cmap='gray_r', vmin=0, vmax=1, interpolation='nearest', extent=[0,1,0,1])

#    red_overlay = ax_left.scatter([], [], color='red', alpha=0, s=30)  # increased size for visibility    
    red_overlay = ax_left.scatter([], [], color='red', alpha=0.4, s=30)  # increased size for visibility
    
    # Epoch boundary rectangle (will update per frame)
    epoch_rect = Rectangle((0,0), 1, 1, linewidth=2, edgecolor='lime', facecolor='none', linestyle='--')
    ax_left.add_patch(epoch_rect)
    
    # Epoch boundary rectangle
#    epoch_rect = Rectangle((0,0), 1, 1, linewidth=2, edgecolor='lime', facecolor='none', linestyle='--')
#    ax_left.add_patch(epoch_rect)
    
    trace_lines = []
    
    def animate(frame):
        amplitudes, N_h, traces = all_data[frame]
        side = math.ceil(math.sqrt(N_h))
        grid = np.zeros((side, side))
        grid.flat[:N_h] = amplitudes / (amplitudes.max() + 1e-8) if amplitudes.max() > 0 else grid
        
        im.set_array(grid)
        im.set_extent([0, side, 0, side])
        ax_left.set_xlim(0, side)
        ax_left.set_ylim(0, side)
        
        # Update epoch boundary
        epoch_rect.set_width(side)
        epoch_rect.set_height(side)
        
        # Clear traces
        for line in trace_lines:
            line.remove()
        trace_lines.clear()
        
        # Red entanglement overlay — now CENTERED in each pixel
        if show_red_intensity:
            marked = np.where(amplitudes > 0)[0]
            valid_marked = marked[marked < N_h]  # strictly within epoch
            if len(valid_marked) > 0:
                amps_marked = amplitudes[valid_marked]
                row_marked = valid_marked // side + 0.5   # ← Center offset!
                col_marked = valid_marked % side + 0.5     # ← Center offset!
                red_overlay.set_offsets(np.column_stack((col_marked, row_marked)))
                # Intensity: higher amplitude → stronger red (darker + larger)
                red_alpha = 0.4 + 0.6 * (amps_marked / amps_marked.max())
                red_overlay.set_alpha(red_alpha.mean())
                red_overlay.set_sizes(20 + 30 * (amps_marked / amps_marked.max()))  # size scales with amplitude
            else:
                red_overlay.set_offsets(np.empty((0,2)))
        else:
            red_overlay.set_offsets(np.empty((0,2)))
        
        if show_traces_panel:
            ax_right.clear()
            ax_right.set_xlim(0, side)
            ax_right.set_ylim(0, side)
            ax_right.set_title("Operator Traces (Crystal Faces)")
            ax_right.set_xlabel("Column")
            ax_right.set_ylabel("Row")
            ax_right.grid(True, alpha=0.2)
        
        # Draw traces
        for x, y0, p, face_z in traces:
            face_idx = R.index(face_z)
            color = face_colors[face_idx]
            t = np.linspace(0, min(10*p, side*1.5), 200)
            wave_y = y0 + t + 0.08 * side * np.sin(2 * np.pi * t / p)
            row = wave_y // side
            col = wave_y % side
            valid = (row >= 0) & (row < side) & (col >= 0) & (col < side)
            if show_traces_panel:
                line, = ax_right.plot(col[valid], row[valid], color=color, alpha=0.5, linewidth=1.2)
            else:
                line, = ax_left.plot(col[valid], row[valid], color=color, alpha=0.1, linewidth=0.6)
            trace_lines.append(line)

        title_left.set_text(f"Epoch {frame+1} | N={N_h} | Holes={np.sum(amplitudes == 0)}")
        artists = [im, red_overlay, epoch_rect] + trace_lines
        return artists

        
       # title_left.set_text(f"Epoch {frame+1} | N={N_h} | Holes={np.sum(amplitudes == 0)}")
       # artists = [im, red_overlay, epoch_rect] + trace_lines
       # return artists

    anim = animation.FuncAnimation(fig, animate, frames=num_epochs, interval=800, blit=False)
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Animation saved to {output_file}")
    print("Red markers are now centered in each pixel square.")

# Run parameters
if __name__ == '__main__':
    classes = [11]  # or [11,13]
    num_epochs = 30
    show_traces_panel = True
    show_red_intensity = False
    
    create_dual_view_animation(classes, num_epochs, show_traces_panel=show_traces_panel, show_red_intensity=show_red_intensity)


"""
def create_deinterlaced_animation(classes, num_epochs, deinterlace_mode=False, output_file='deinterlaced_epoch_animation.gif'):
    pool = mp.Pool(mp.cpu_count())
    args_list = [(classes, h) for h in range(1, num_epochs + 1)]
    all_data = pool.map(precompute_epoch, args_list)
    pool.close()
    
    fig, ax = plt.subplots(figsize=(12, 12))
    title = ax.set_title("Epoch 1 – Amplitude Field")
    im = ax.imshow(np.zeros((1,1)), cmap='gray_r', vmin=0, vmax=1, interpolation='nearest', extent=[0,1,0,1])
    wave_lines = []
    red_overlay = ax.scatter([], [], color='red', alpha=0.3, s=20)
    
    if deinterlace_mode:
        # Deinterlace into 24 residue classes (one per sub-animation frame cycle)
        num_frames = num_epochs * 24
    else:
        num_frames = num_epochs
    
    def animate(frame):
        if deinterlace_mode:
            epoch_idx = frame // 24
            residue_idx = frame % 24
            amplitudes, N_h, traces = all_data[epoch_idx]
            # Filter traces to only this residue class
            residue_z = R[residue_idx]
            filtered_traces = [t for t in traces if t[3] == residue_z]
            title_suffix = f"Residue Class {residue_z}"
        else:
            epoch_idx = frame
            amplitudes, N_h, traces = all_data[epoch_idx]
            filtered_traces = traces
            title_suffix = ""
        
        side = math.ceil(math.sqrt(N_h))
        grid = np.zeros((side, side))
        grid.flat[:N_h] = amplitudes / (amplitudes.max() + 1e-8) if amplitudes.max() > 0 else grid
        
        im.set_array(grid)
        im.set_extent([0, side, 0, side])
        ax.set_xlim(0, side)
        ax.set_ylim(0, side)
        
        # Clear waves
        for line in wave_lines:
            line.remove()
        wave_lines.clear()
        
        # Draw filtered wave fronts for this deinterlaced class
        for x, y0, p, face_z in filtered_traces:
            face_idx = R.index(face_z)
            color = face_colors[face_idx]
            t = np.linspace(0, min(10*p, side*1.5), 200)
            wave_y = y0 + t + 0.08 * side * np.sin(2 * np.pi * t / p)
            row = wave_y // side
            col = wave_y % side
            valid = (row >= 0) & (row < side) & (col >= 0) & (col < side)
            line, = ax.plot(col[valid], row[valid], color=color, alpha=0.5, linewidth=1.2)
            wave_lines.append(line)
        
        # Red overlay (global for amplitude, not filtered)
        marked = np.where(amplitudes > 0)[0]
        if len(marked) > 0:
            amps_marked = amplitudes[marked]
            row_marked = marked // side
            col_marked = marked % side
            red_overlay.set_offsets(np.column_stack((col_marked, row_marked)))
            red_alpha = 0.3 + 0.7 * (amps_marked / amps_marked.max())
            red_overlay.set_alpha(red_alpha.mean())
        
        title.set_text(f"Epoch {epoch_idx+1} | N={N_h} | Holes={np.sum(amplitudes == 0)} | {title_suffix}")
        return [im, red_overlay] + wave_lines

    anim = animation.FuncAnimation(fig, animate, frames=num_frames, interval=500, blit=True)
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Animation saved to {output_file}")
    print(f"Deinterlace mode is {'ON' if deinterlace_mode else 'OFF'}")

# Run parameters
if __name__ == '__main__':
    classes = [11]  # or [11,13]
    num_epochs = 30
    deinterlace_mode = True  # <--- Toggle here (True = deinterlace into 24 residue sequences)
    create_deinterlaced_animation(classes, num_epochs, deinterlace_mode=deinterlace_mode)
    """