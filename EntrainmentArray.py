import math
import cmath
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

def epoch_length(h):
    return 90 * h**2 - 12 * h + 1

def simulate_entrainment_diagonals(h_max=5, grid_size=20, base_start=0):
    class11 = [11]
    class13 = [13]
    results = []
    
    current_start = base_start
    for h in range(1, h_max + 1):
        epoch_len = epoch_length(h)
        amps11 = compute_amplitudes(class11, current_start, epoch_len)
        amps13 = compute_amplitudes(class13, current_start, epoch_len)
        
        # For visualization, take first grid_size**2 elements, reshape to grid
        total = min(epoch_len, grid_size * grid_size)
        grid11 = np.zeros((grid_size, grid_size))
        grid13 = np.zeros((grid_size, grid_size))
        grid_joint = np.zeros((grid_size, grid_size))
        
        for i in range(total):
            r, c = divmod(i, grid_size)
            marked11 = amps11[i] > 0
            marked13 = amps13[i] > 0
            grid11[r, c] = 1 if marked11 else 0
            grid13[r, c] = 1 if marked13 else 0
            grid_joint[r, c] = 1 if marked11 or marked13 else 0  # 1 if composite in at least one
        
        # Save plots
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))
        axs[0].imshow(grid11, cmap='binary')
        axs[0].set_title(f'Class 11 Epoch {h}')
        axs[1].imshow(grid13, cmap='binary')
        axs[1].set_title(f'Class 13 Epoch {h}')
        axs[2].imshow(grid_joint, cmap='binary')
        axs[2].set_title(f'Joint Epoch {h}')
        plt.savefig(f'epoch_{h}.png')
        plt.close()
        
        # Text representation for output
        text11 = '\n'.join(' '.join('X' if cell else '.' for cell in row) for row in grid11)
        text13 = '\n'.join(' '.join('X' if cell else '.' for cell in row) for row in grid13)
        text_joint = '\n'.join(' '.join('X' if cell else '.' for cell in row) for row in grid_joint)
        
        results.append(f"Epoch {h} (length {epoch_len}, start {current_start})\nClass 11:\n{text11}\n\nClass 13:\n{text13}\n\nJoint:\n{text_joint}\n")
        
        current_start += epoch_len
    
    return '\n'.join(results)

# Run simulation
output = simulate_entrainment_diagonals(h_max=3, grid_size=10, base_start=0)
print(output)