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

def compute_epoch_marks(classes, x_max, min_index=0):
    """
    Compute all marked indices up to epoch x_max.
    Returns a list of all unique marked indices (sorted) across epochs 1 to x_max.
    """
    all_ops = []
    for k in classes:
        all_ops.extend(get_operators(k))
    
    marked_set = set()
    
    for x in range(1, x_max + 1):
        for l_val, m_val, z, primitive in all_ops:
            y0 = 90 * x * x - l_val * x + m_val
            p = z + 90 * (x - 1)
            if p <= 0:
                continue
            current = y0
            while current >= min_index:
                marked_set.add(current)
                current += p
                if current > 10**18:  # Prevent infinite loop in large x
                    break
    
    return sorted(list(marked_set))

def create_epoch_animation(classes, max_epochs, output_file='epoch_crystal_growth.gif'):
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(np.zeros((10, 10, 3)), interpolation='nearest')  # Placeholder
    title_text = ax.set_title("Epoch 0: Initial (empty)")
    ax.set_xlabel("Column")
    ax.set_ylabel("Row")
    ax.axis('off')
    
    def animate(epoch):
        # Compute cumulative marked indices up to this epoch
        marked_indices = compute_epoch_marks(classes, epoch)
        
        if not marked_indices:
            grid_rgb = np.zeros((10, 10, 3))
            ax.set_title(f"Epoch {epoch}: No marks yet")
            im.set_array(grid_rgb)
            return [im]
        
        # Find the smallest perfect square >= max index
        max_idx = max(marked_indices)
        n = math.ceil(math.sqrt(max_idx + 1))  # +1 for safety
        total_pixels = n * n
        
        # Create amplitude grid (0 = hole, 1 = marked)
        amp_grid = np.zeros(total_pixels)
        for idx in marked_indices:
            if idx < total_pixels:
                amp_grid[idx] = 1  # Simple binary for now: marked or not
        
        # Reshape to n x n
        grid = amp_grid.reshape(n, n)
        
        # Create RGB: black holes, bright white marks
        grid_rgb = np.stack([grid, grid, grid], axis=-1)  # Grayscale
        # Highlight holes (0) in yellow
        holes = (grid == 0)
        hole_rgb = plt.cm.colors.to_rgba('#FFFF00')[:3]
        for c in range(3):
            grid_rgb[:,:,c][holes] = hole_rgb[c]
        
        im.set_array(grid_rgb)
        ax.set_title(f"Epoch {epoch} - Grid {n}×{n} (Total indices covered: ~{n*n})")
        ax.set_xlim(0, n)
        ax.set_ylim(0, n)
        
        return [im]
    
    anim = animation.FuncAnimation(fig, animate, frames=range(1, max_epochs + 1),
                                   interval=8, blit=True)
    
    anim.save(output_file, writer='pillow', fps=1)
    plt.close(fig)
    
    print(f"Epoch-growing animation saved to {output_file}")

# ================= CONFIGURABLE PARAMETERS =================
classes = [11]      # Example class
max_epochs = 1     # Number of epochs to animate (adjust for speed vs detail)

# ================= RUN =================
if __name__ == '__main__':
    create_epoch_animation(classes, max_epochs)