import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

def get_operators(k):
    R = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89]
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

def compute_amplitudes(classes, start, width):
    segment_start = start
    segment_end = start + width
    a = 90
    b = -300
    c = 250 - segment_end
    d = (b**2) - (4 * a * c)
    sol2 = (-b + math.sqrt(d)) / (2 * a) if d >= 0 else (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = int(sol2.real) + 1
    
    list_amp = [0] * (segment_end - segment_start)
    
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
    
    return list_amp

def generate_pixel_map(amplitudes, n_rows, m_cols, output_file='amplitude_map.png'):
    total_pixels = n_rows * m_cols
    # Pad or truncate amplitudes to fit grid
    amps = amplitudes[:total_pixels]
    if len(amps) < total_pixels:
        amps += [0] * (total_pixels - len(amps))  # Pad with black (holes)
    
    # Normalize to [0, 255] for grayscale brightness
    max_amp = max(amps) if max(amps) > 0 else 1
    normalized = [int(255 * (amp / max_amp)) for amp in amps]
    
    # Reshape into n x m grid
    grid = np.array(normalized).reshape(n_rows, m_cols)
    
    # Plot as grayscale image (0=black, 255=white/bright)
    plt.figure(figsize=(m_cols/100, n_rows/100))  # Adjust size for visibility
    plt.imshow(grid, cmap='gray', vmin=0, vmax=255)
    plt.axis('off')  # No axes for clean pixel map
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0)
    plt.close()
    print(f"Pixel map saved to {output_file}")

# Configurable parameters
classes = [11]  # Example: class 11
start_index = 0  # Starting index
segment_width = 81000000  # Total amplitudes to compute (must be n_rows * m_cols)
n_rows = 9000  # Number of rows in grid
m_cols = 9000  # Number of columns in grid (width = n_rows * m_cols = 10000)

# Compute amplitudes
amplitudes = compute_amplitudes(classes, start_index, segment_width)

# Generate and save the pixel map
generate_pixel_map(amplitudes, n_rows, m_cols)