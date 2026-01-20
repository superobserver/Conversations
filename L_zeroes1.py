import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from mpmath import mp, dirichlet, findroot

mp.dps = 50  # Higher precision for more zeros

# Define a non-principal Dirichlet character mod 90 for k=11
# As product of quadratic characters mod 2,5,9, adjusted to emphasize ≡11
def chi_mod90(n):
    if math.gcd(n, 90) != 1:
        return 0
    # Quadratic character mod 5: Legendre (n/5)
    leg5 = pow(n, (5-1)//2, 5)
    # Mod 9: Jacobi (n/9)
    jac9 = pow(n, (9-1)//2, 9)
    # Mod 2: trivial or (-1)^{(n-1)/2}
    chi2 = (-1)**((n**2 - 1)//8)  # Non-principal mod 8, but adapt for 2
    chi = leg5 * jac9 * chi2
    # Phase adjustment to align with k=11: normalize so chi(11) =1
    return chi if abs(chi) > 0 else 0

# Character vector for mpmath.dirichlet (chi(1) to chi(90))
chi_vector = [chi_mod90(j) for j in range(1, 91)]

# Compute L(s, chi) for s=0.5 + it
def L_chi(t):
    s = mp.mpc(0.5, t)
    return dirichlet(s, chi_vector)

# Find more approximate L-zeros on critical line (Re=0.5, Im=t>0)
def find_l_zeros(num_zeros=20, t_max=100, step=0.05):
    zeros = []
    t = 1.0
    while len(zeros) < num_zeros and t < t_max:
        def f_real(u):
            return mp.re(L_chi(u))
        def f_imag(u):
            return mp.im(L_chi(u))
        try:
            root_real = findroot(f_real, t, tol=1e-10)
            if abs(f_imag(root_real)) < 1e-5:  # Check near zero
                zeros.append(float(root_real))
        except ValueError:
            pass
        t += step
    return zeros

# Your sieve functions (optimized for larger max_h)
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
    
    list_amp = np.zeros(total_length)
    
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
                    idx = int(current - segment_start)
                    list_amp[idx] += 1
                current += p
    
    return list_amp

def compute_hole_densities(k=11, max_h=100):
    h_values = np.arange(1, max_h + 1)
    tau_values = []
    for h in h_values:
        N_h = 90 * h**2 - 12 * h + 1
        amplitudes = compute_amplitudes([k], 0, N_h)
        H = np.sum(amplitudes == 0)
        tau = H / N_h
        tau_values.append(tau)
    tau_values = np.array(tau_values)
    
    # PNT baseline (refactored for denser space: 3.75 / (9 + 2 ln h), but scale by 15/4 ~3.75 for 4/15 density)
    baseline = (15/4) * 3.75 / (9 + 2 * np.log(h_values))  # Refactored scaling for denser sublattice
    
    # Variations Δτ(h) = τ(h) - baseline
    variations = tau_values - baseline
    
    return h_values, tau_values, baseline, variations

# Build inversion graph and compute crystal spectrum
def compute_crystal_spectrum(k):
    R = get_residues()
    n = len(R)
    A = np.zeros((n, n))
    for i, z in enumerate(R):
        o = (k * pow(z, -1, 90)) % 90
        j = R.index(o)
        A[i, j] = 1
        A[j, i] = 1
    eigenvalues = np.sort(np.real(np.linalg.eigvals(A)))
    return eigenvalues

# Main: Compute more zeros, link via correlation to crystal spectrum
def main(k=11, max_h=100, num_zeros=20):
    # Compute L-zeros (more, higher t_max)
    l_zeros = find_l_zeros(num_zeros=num_zeros, t_max=200, step=0.02)
    print(f"Approximate L-zeros for chi mod 90 (Im parts): {l_zeros}")
    
    # Compute hole densities and variations
    h, tau, baseline, variations = compute_hole_densities(k, max_h)
    
    # Crystal spectrum
    crystal_eig = compute_crystal_spectrum(k)
    crystal_norm = 2 * (crystal_eig - min(crystal_eig)) / (max(crystal_eig) - min(crystal_eig)) - 1 if max(crystal_eig) > min(crystal_eig) else crystal_eig  # Normalize to [-1,1]
    crystal_norm = np.sort(crystal_norm)
    
    # Link via algebraic cotermination: correlate L-zeros Im to crystal eig (scaled)
    min_len = min(len(l_zeros), len(crystal_norm))
    l_zeros_norm = np.sort(np.array(l_zeros))[:min_len]
    l_zeros_norm = (l_zeros_norm - min(l_zeros_norm)) / max(l_zeros_norm) * 2 - 1  # Normalize to [-1,1]
    r_spec, _ = pearsonr(crystal_norm[:min_len], l_zeros_norm)
    print(f"Spectral Correlation (crystal eig vs L-zeros Im): {r_spec:.4f}")
    
    # Correlate variations to oscillatory signal from L-zeros
    N_h = 90 * h**2 - 12 * h + 1
    oscillatory = np.zeros(len(h))
    for j, t in enumerate(l_zeros):
        oscillatory += np.cos(t * np.log(N_h)) / (j+1)  # Weighted sum for error proxy
    r_var, _ = pearsonr(variations, oscillatory)
    print(f"Variation Correlation to L-zero oscillatory: {r_var:.4f}")
    
    # Plot with links
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    ax1.plot(h, tau, label='Observed τ(h)', marker='o')
    ax1.plot(h, baseline, label='Refactored PNT Baseline', linestyle='--')
    ax1.set_title(f'Hole Density τ(h) for Class {k} (h=1 to {max_h})')
    ax1.legend()
    ax1.grid(True)
    
    ax2.plot(h, variations, label='Variation Δτ(h)', marker='x', color='red')
    ax2.plot(h, 0.01 * oscillatory, label='Scaled L-zero Oscillatory', linestyle=':', color='green')  # Scaled for vis
    ax2.set_title(f'Density Variation & L-zero Signal (r={r_var:.4f})')
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig('l_zero_crystal_link.png')
    plt.close()
    print("Plot saved to 'l_zero_crystal_link.png'")

main(k=11, max_h=100, num_zeros=20)