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

def compute_cumulative_marks(classes, start_index, max_epoch, truncate_to_square=True):
    """
    Compute cumulative marked indices up to max_epoch.
    Returns list of lists: marks_per_epoch[x-1] = list of marked indices in epoch x only.
    """
    all_ops = []
    for k in classes:
        all_ops.extend(get_operators(k))
    
    marks_per_epoch = []
    cumulative_marks = set()
    
    for x in range(1, max_epoch + 1):
        epoch_marks = set()
        for l_val, m_val, z, primitive in all_ops:
            y0 = 90 * x * x - l_val * x + m_val
            p = z + 90 * (x - 1)
            if p <= 0:
                continue
            current = y0
            while current >= start_index:
                if current not in cumulative_marks:
                    epoch_marks.add(current)
                current += p
                if current > start_index + 10**12:  # Safety bound
                    break
        marks_per_epoch.append(sorted(epoch_marks))
        cumulative_marks.update(epoch_marks)
    
    return marks_per_epoch

def create_epoch_growing_animation(classes, start_index, max_epochs, output_file='epoch_growing_crystal.gif'):
    marks_per_epoch = compute_cumulative_marks(classes, start_index, max_epochs)
    
    fig, ax = plt.subplots(figsize=(10, 10))
    # Initial empty placeholder
    initial_n = 10
    im = ax.imshow(np.zeros((initial_n, initial_n, 3)), interpolation='nearest')
    title_text = ax.set_title("Epoch 0: Initial (empty)")
    ax.set_xlabel("Column")
    ax.set_ylabel("Row")
    ax.axis('off')
    
    def animate(epoch):
        # Cumulative marks up to this epoch
        all_marks = set()
        for e in range(epoch):
            all_marks.update(marks_per_epoch[e])
        
        if not all_marks:
            grid_rgb = np.zeros((10, 10, 3))
            ax.set_title(f"Epoch {epoch}: No marks yet")
            im.set_array(grid_rgb)
            return [im]
        
        max_idx = max(all_marks)
        n = math.ceil(math.sqrt(max_idx + 1))
        if n % 2 != 0:
            n += 1  # Ensure even for perfect square
        
        total_pixels = n * n
        
        # Binary grid: 1 = marked, 0 = hole
        amp_grid = np.zeros(total_pixels)
        for idx in all_marks:
            if idx < total_pixels:
                amp_grid[idx] = 1
        
        grid = amp_grid.reshape(n, n)
        
        # RGB: black holes, bright white marks
        grid_rgb = np.stack([grid, grid, grid], axis=-1)
        holes = (grid == 0)
        hole_rgb = plt.cm.colors.to_rgba('#FFFF00')[:3]
        for c in range(3):
            grid_rgb[:,:,c][holes] = hole_rgb[c]
        
        im.set_array(grid_rgb)
        ax.set_xlim(0, n)
        ax.set_ylim(0, n)
        ax.set_title(f"Epoch {epoch} - Grid {n}×{n} (Cumulative marks: {len(all_marks)})")
        return [im]
    
    anim = animation.FuncAnimation(fig, animate, frames=range(max_epochs + 1),
                                   interval=8, blit=True)
    
    anim.save(output_file, writer='pillow', fps=1)
    plt.close(fig)
    
    print(f"Epoch-growing animation saved to {output_file}")

# ================= CONFIGURABLE PARAMETERS =================
classes = [11]      # Example class
start_index = 0     # Starting index
max_epochs = 3     # Number of epochs to animate (adjust for detail vs speed)

# ================= RUN =================
if __name__ == '__main__':
    create_epoch_growing_animation(classes, start_index, max_epochs)