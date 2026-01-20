import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
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

def create_cylinder_animation(classes, num_epochs, cylinder_mode=True, output_file='cylinder_epoch_animation.gif'):
    pool = mp.Pool(mp.cpu_count())
    args_list = [(classes, h) for h in range(1, num_epochs + 1)]
    all_data = pool.map(precompute_epoch, args_list)
    pool.close()
    
    fig = plt.figure(figsize=(12, 12))
    if cylinder_mode:
        ax = fig.add_subplot(111, projection='3d')
    else:
        ax = fig.add_subplot(111)
    
    title = fig.suptitle("Epoch 1 – Amplitude Field")
    
    def animate(frame):
        amplitudes, N_h, traces = all_data[frame]
        side = math.ceil(math.sqrt(N_h))
        grid = np.zeros((side, side))
        grid.flat[:N_h] = amplitudes / (amplitudes.max() + 1e-8) if amplitudes.max() > 0 else grid
        
        ax.clear()
        if cylinder_mode:
            # Cylinder mapping: x (col) → θ, y (row) → z height, r=1
            theta = 2 * np.pi * np.arange(side) / side
            z_height = np.arange(side)
            Theta, Z = np.meshgrid(theta, z_height)
            R = np.ones_like(Theta)
            X = R * np.cos(Theta)
            Y = R * np.sin(Theta)
            C = grid.T  # Transpose for correct orientation
            ax.plot_surface(X, Y, Z, facecolors=plt.cm.gray_r(C), rstride=1, cstride=1, linewidth=0, alpha=0.8)
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z (Height)')
            ax.set_title('Cylinder Mapping')
            # Rotate view for dynamic effect
            ax.view_init(elev=20, azim=frame * 4)  # Gentle rotation
        else:
            im = ax.imshow(grid, cmap='gray_r', vmin=0, vmax=1, interpolation='nearest', extent=[0, side, 0, side])
            ax.set_xlim(0, side)
            ax.set_ylim(0, side)
            ax.set_title('Square Mapping')
        
        title.set_text(f"Epoch {frame+1} | N={N_h} | Holes={np.sum(amplitudes == 0)}")
        return []

    anim = animation.FuncAnimation(fig, animate, frames=num_epochs, interval=500, blit=False)
    anim.save(output_file, writer='pillow')
    plt.close(fig)
    
    print(f"Animation saved to {output_file}")

# Run parameters
if __name__ == '__main__':
    classes = [11]  # Example class
    num_epochs = 40  # Number of epochs
    cylinder_mode = True  # Set to False for square mapping
    create_cylinder_animation(classes, num_epochs, cylinder_mode=cylinder_mode)