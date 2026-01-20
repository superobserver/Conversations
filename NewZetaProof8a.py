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
    red_overlay = ax_left.scatter([], [], color='red', alpha=0, s=30)  # increased size for visibility
    
    # Epoch boundary rectangle
    epoch_rect = Rectangle((0,0), 1, 1, linewidth=2, edgecolor='lime', facecolor='none', linestyle='--')
    ax_left.add_patch(epoch_rect)
    
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
            valid_marked = marked[marked < N_h]
            if len(valid_marked) > 0:
                amps_marked = amplitudes[valid_marked]
                row_marked = valid_marked // side + 0.5   # ← Center offset!
                col_marked = valid_marked % side + 0.5     # ← Center offset!
                red_overlay.set_offsets(np.column_stack((col_marked, row_marked)))
                # Intensity: higher amplitude → stronger red
                red_alpha = 0.3 + 0.7 * (amps_marked / amps_marked.max())
                red_overlay.set_alpha(red_alpha.mean())
                # Size also scales with amplitude for emphasis
                red_overlay.set_sizes(20 + 30 * (amps_marked / amps_marked.max()))
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
    show_red_intensity = True
    
    create_dual_view_animation(classes, num_epochs, show_traces_panel=show_traces_panel, show_red_intensity=show_red_intensity)