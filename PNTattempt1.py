import math
import cmath
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import multiprocessing as mp

import math
import cmath
import numpy as np

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

def precompute_epoch(class_k, h):
    N_h = 90 * h**2 - 12 * h + 1
    amplitudes = compute_amplitudes([class_k], 0, N_h)
    hole_count = np.sum(amplitudes == 0)
    return hole_count

def plot_hole_oscillations(class_k=11, max_h=1000, sma_window=50):
    h_values = np.arange(1, max_h + 1)
    hole_counts = []
    
    with mp.Pool(mp.cpu_count()) as pool:
        hole_counts = pool.starmap(precompute_epoch, [(class_k, h) for h in h_values])
    
    hole_counts = np.array(hole_counts)
    
    # Moving average (SMA)
    moving_avg = np.convolve(hole_counts, np.ones(sma_window)/sma_window, mode='valid')
    padded_moving_avg = np.pad(moving_avg, (sma_window-1, 0), mode='constant', constant_values=np.nan)
    
    # Undulation delta = H - avg
    delta = hole_counts - padded_moving_avg
    
    # Rate of change of undulation (first difference)
    rate_change = np.diff(delta)
    rate_change = np.insert(rate_change, 0, np.nan)  # Align lengths
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    ax1.plot(h_values, hole_counts, label='Hole Count \( H_k(h) \)', color='blue')
    ax1.plot(h_values, padded_moving_avg, label=f'SMA (window={sma_window})', color='orange', linestyle='--')
    ax1.set_xlabel('Epoch \( h \)')
    ax1.set_ylabel('Hole Count')
    ax1.set_title(f'Hole Count vs. Epoch for Class {class_k} (h=1 to {max_h})')
    ax1.legend()
    ax1.grid(True)
    
    ax2.plot(h_values, rate_change, label='Rate of Change of Undulation \( \\Delta\'(h) \)', color='green')
    ax2.set_xlabel('Epoch \( h \)')
    ax2.set_ylabel('Rate of Change')
    ax2.set_title('Signal of Undulation Rate')
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig('hole_oscillation_analysis.png')
    plt.close(fig)
    
    print(f"Analysis complete. Plot saved as 'hole_oscillation_analysis.png'. Observed crossings: {np.sum(np.diff(np.sign(delta[1:])) != 0)}")

# Run the extension


plot_hole_oscillations(class_k=11, max_h=1000)