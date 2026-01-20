import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# First 24 imaginary parts of non-trivial zeta zeros (from Odlyzko et al.)
ZETA_ZERO_IMAG = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178, 40.918719, 43.327073,
                  48.005151, 49.773832, 52.970321, 56.446248, 59.347045, 60.831778, 65.112544, 67.079811,
                  69.546402, 72.067158, 75.704690, 79.337375, 82.910381, 84.735493, 87.425274, 88.809111]

def get_residues():
    return [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]

def build_inversion_graph(k):
    R = get_residues()
    n = len(R)
    A = np.zeros((n, n))
    for i, z in enumerate(R):
        o = (k * pow(z, -1, 90)) % 90
        j = R.index(o)
        A[i, j] = 1
        A[j, i] = 1  # Undirected
    return A

def compute_crystal_spectrum(k):
    A = build_inversion_graph(k)
    eigenvalues = np.sort(np.real(np.linalg.eigvals(A)))
    return eigenvalues

def get_zero_periods(scale=10):
    # Scale and round to integers for period use
    return np.round(np.array(ZETA_ZERO_IMAG) * scale).astype(int)

def get_operators(k, zero_interlace=False):
    R = get_residues()
    if zero_interlace:
        R = get_zero_periods() % 90  # Mod 90 to fit residue space
        R = [r if r in R else R[r % 24] for r in R]  # Map to valid residues if needed
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

def compute_amplitudes(classes, h, zero_interlace=False):
    N_h = 90 * h**2 - 12 * h + 1
    list_amp = np.zeros(N_h)
    all_ops = []
    for k in classes:
        all_ops.extend(get_operators(k, zero_interlace))
    a = 90
    b = -300
    c = 250 - N_h
    d = (b**2) - (4 * a * c)
    sol2 = (-b + math.sqrt(d)) / (2 * a) if d >= 0 else (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = int(sol2.real) + 1
    for x in range(1, new_limit + 1):
        for l_val, m_val, z, primitive in all_ops:
            y0 = 90 * x * x - l_val * x + m_val
            p = z + 90 * (x - 1)
            if p <= 0 or y0 >= N_h:
                continue
            current = max(y0, 0)
            while current < N_h:
                list_amp[int(current)] += 1
                current += p
    return list_amp

def compute_hole_count(amplitudes):
    return np.sum(amplitudes == 0)

def compare_spectra(k=11, h=10):
    # Crystal spectrum
    crystal_eig = compute_crystal_spectrum(k)
    
    # Zeta zero imaginaries (normalized to [-1,1] range for comparison)
    zeta_imag = np.array(ZETA_ZERO_IMAG)
    zeta_norm = 2 * (zeta_imag - min(zeta_imag)) / (max(zeta_imag) - min(zeta_imag)) - 1  # to [-1,1]
    zeta_norm = np.sort(zeta_norm)
    
    # Correlation
    min_len = min(len(crystal_eig), len(zeta_norm))
    correlation, _ = pearsonr(crystal_eig[:min_len], zeta_norm[:min_len])
    
    # Plot comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(range(min_len), crystal_eig[:min_len], label='Crystal Eigenvalues', marker='o')
    ax.plot(range(min_len), zeta_norm[:min_len], label='Normalized Zeta Zeros Im', marker='x')
    ax.set_title(f'Spectra Comparison for Class {k} (Correlation: {correlation:.4f})')
    ax.set_xlabel('Index')
    ax.set_ylabel('Value')
    ax.legend()
    ax.grid(True)
    plt.savefig('spectra_comparison.png')
    plt.close()

    # Sieve comparison
    amps_original = compute_amplitudes([k], h)
    amps_zero = compute_amplitudes([k], h, zero_interlace=True)
    
    holes_original = compute_hole_count(amps_original)
    holes_zero = compute_hole_count(amps_zero)
    
    N_h = len(amps_original)
    print(f"Epoch h={h}, N={N_h}")
    print(f"Original Sieve: Holes={holes_original}, τ={holes_original/N_h:.4f}")
    print(f"Zero-Interlaced Sieve: Holes={holes_zero}, τ={holes_zero/N_h:.4f}")
    print("Spectra plot saved to 'spectra_comparison.png'")

# Run the comparison
compare_spectra(k=11, h=10)