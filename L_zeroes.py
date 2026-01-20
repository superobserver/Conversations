import math
import cmath
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from mpmath import mp, dirichlet, findroot

mp.dps = 30  # High precision for L-function evaluation

# Define a non-principal Dirichlet character mod 90 for k=11
# Mod 90 = 2*5*9, characters product of mod2, mod5, mod9
# Example: chi(n) = (n/5)_legendre * chi_9(n), adjusted to emphasize ≡11 mod 90
def dirichlet_char_mod90(n):
    if math.gcd(n, 90) != 1:
        return 0
    # Legendre symbol mod 5 (quadratic character)
    leg5 = pow(n, (5-1)//2, 5)  # (n/5)
    # Mod 9: primitive character, e.g., chi(1)=1, chi(2)=i, but simplify to real
    # For simplicity, use a real non-principal: chi(n) = (n/5) * (n/9) where (·/9) is Jacobi
    jac9 = pow(n, (9-1)//2, 9) if 9 > 2 else 1
    chi = leg5 * jac9
    # Adjust phase to align with k=11: if n≡11 mod 90, chi should be non-zero
    return chi if chi else 0  # Filter non-coprime

# Compute L(s, chi) for s=1/2 + it
def L_chi(t):
    s = mp.mpc(0.5, t)
    # mpmath.dirichletl(s, chi, q=90) but chi as list on generators
    # Generators for mod 90: e.g., 11, 17 (primitive?)
    # For simplicity, approximate with mp.dirichlet(s, chi_vector) but define chi_vector as list chi(1) to chi(90)
    chi_vector = [dirichlet_char_mod90(j) for j in range(1, 91)]
    return dirichlet(s, chi_vector)

# Find approximate L-zeros on critical line (Re=0.5, Im=t>0) by root-finding
def find_l_zeros(num_zeros=8, t_max=50, step=0.1):
    zeros = []
    t = 1.0
    while len(zeros) < num_zeros and t < t_max:
        # Find where Re(L) crosses zero (phase change or numerical)
        def f_real(u):
            return mp.re(L_chi(u))
        try:
            root = findroot(f_real, t)
            if abs(mp.im(L_chi(root))) < 1e-5:  # Check near zero
                zeros.append(float(root))
        except ValueError:
            pass
        t += step
    return zeros

# Your sieve functions (from previous)
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

# Compute hole density τ(h) and variations for h=1 to max_h
def compute_hole_densities(k=11, max_h=50):
    h_values = np.arange(1, max_h + 1)
    tau_values = []
    for h in h_values:
        N_h = 90 * h**2 - 12 * h + 1
        amplitudes = compute_amplitudes([k], 0, N_h)
        H = np.sum(amplitudes == 0)
        tau = H / N_h
        tau_values.append(tau)
    tau_values = np.array(tau_values)
    
    # PNT baseline (expected τ ≈ 3.75 / (9 + 2 ln h))
    baseline = 3.75 / (9 + 2 * np.log(h_values))
    
    # Variations Δτ(h) = τ(h) - baseline
    variations = tau_values - baseline
    
    return h_values, tau_values, baseline, variations

# Main computation and correlation
def main(k=11, max_h=50, num_zeros=8):
    # Compute L-zeros
    l_zeros = find_l_zeros(num_zeros=num_zeros)
    print(f"Approximate L-zeros for chi mod 90 (Im parts): {l_zeros}")
    
    # Compute hole densities and variations
    h, tau, baseline, variations = compute_hole_densities(k, max_h)
    
    # Plot density and variation
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    ax1.plot(h, tau, label='Observed τ(h)', marker='o')
    ax1.plot(h, baseline, label='PNT Baseline', linestyle='--')
    ax1.set_title(f'Hole Density τ(h) for Class {k} (h=1 to {max_h})')
    ax1.set_xlabel('Epoch h')
    ax1.set_ylabel('Density τ')
    ax1.legend()
    ax1.grid(True)
    
    ax2.plot(h, variations, label='Variation Δτ(h)', marker='x', color='red')
    ax2.set_title('Density Variation Δτ(h)')
    ax2.set_xlabel('Epoch h')
    ax2.set_ylabel('Δτ')
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig('hole_density_correlation.png')
    plt.close()
    
    # Correlate variations to cos(t ln N_h) for each zero t
    N_h = 90 * h**2 - 12 * h + 1
    correlations = []
    for t in l_zeros:
        oscillatory = np.cos(t * np.log(N_h))
        r, _ = pearsonr(variations, oscillatory)
        correlations.append(r)
        print(f"Correlation with cos({t} ln N_h): {r:.4f}")
    
    print("Plot saved to 'hole_density_correlation.png'")

main(k=11, max_h=50, num_zeros=8)