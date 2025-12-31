import numpy as np
import matplotlib.pyplot as plt

# Approximate li(x) vectorized (better than sympy for speed)
def li_approx(N):
    mask = N > 1
    result = np.zeros_like(N, dtype=float)
    logN = np.log(N[mask])
    result[mask] = N[mask] / logN + N[mask] / logN**2  # 2-term series
    return result

# Sieve for class 17
h = 50  # Increase as needed (h=3 for quick test)
epoch = 90 * h * h - 12 * h + 1
print(f"Epoch: {epoch}")

params_17 = [
    (132,45,7),(102,20,11),(138,52,13),(72,-1,17),(108,29,19),(78,8,23),
    (138,52,29),(102,28,31),(72,11,37),(132,45,41),(78,16,43),(102,28,47),
    (48,3,49),(108,29,53),(78,16,59),(42,4,61),(102,20,67),(72,11,71),
    (18,0,73),(42,4,77),(78,8,79),(48,3,83),(18,0,89),(72,-1,91)
]

amplitude_lists = [[0] * epoch for _ in range(24)]
a, b, c = 90, -300, 250 - epoch
new_limit = int(((-b + np.sqrt(b**2 - 4*a*c)) / (2*a)).real)

def mark(x, l, m, z, lst):
    y = 90*x*x - l*x + m
    if y >= len(lst): return
    lst[int(y)] += 1
    p = z + 90*(x-1)
    n = 1
    while y + p*n < len(lst):
        lst[int(y + p*n)] += 1
        n += 1

for k, (l,m,z) in enumerate(params_17):
    for x in range(1, new_limit + 1):
        mark(x, l, m, z, amplitude_lists[k])

class17_amp = np.array([sum(amplitude_lists[k][n] for k in range(24)) for n in range(epoch)])

mean_amp = class17_amp.mean()
print(f"Class 17 mean amplitude: {mean_amp:.4f}")

# Exact prime counting
primes17_exact = np.cumsum(class17_amp == 0)

# Main term
n_values = np.arange(epoch)
N_values = 90 * n_values + 17
approx_main = li_approx(N_values) / 24

# Fourier error correction
signal_demean = class17_amp - mean_amp
Y = np.fft.fft(signal_demean)
freq = np.fft.fftfreq(epoch)
K = 50  # Top frequencies
sorted_indices = np.argsort(np.abs(Y))[::-1][:K]
approx_demean = np.zeros(epoch)
for idx in sorted_indices:
    approx_demean += (Y[idx] * np.exp(2j * np.pi * freq[idx] * n_values)).real / epoch

error_approx = -approx_demean / mean_amp
pi_approx = approx_main + np.cumsum(error_approx)

# Plot
plt.figure(figsize=(12,6))
plt.plot(n_values, primes17_exact, label='Exact π_class(n) from Sieve', lw=2)
plt.plot(n_values, approx_main, '--', label='Main Term: li(90n+17)/24', lw=2)
plt.plot(n_values, pi_approx, '-.', label='Approx with Fourier Error Correction', lw=2)
plt.xlabel('n (Index)')
plt.ylabel('Number of Primes in Class 17 mod 90')
plt.title('Prime Counting Function for Class 17 mod 90')
plt.legend()
plt.grid(alpha=0.3)
plt.show()

# Samples
print("Sample: n=1000, Exact:", primes17_exact[min(1000, epoch-1)], "Main Approx:", approx_main[min(1000, epoch-1)], "Fourier Approx:", pi_approx[min(1000, epoch-1)])
print("Sample: n=epoch-1, Exact:", primes17_exact[-1], "Main Approx:", approx_main[-1], "Fourier Approx:", pi_approx[-1])

# ... (previous sieve and FFT setup for class17_amp)
# Assume s_grid and Z17 from continued_modular_zeta as in FourierGrok3.py

# Find pseudo-zeros: minima of |Z17| on critical line
t_values = s_grid.imag
magnitudes = np.abs(Z17)
# Sort by magnitude, take top 10 minima (excluding trivial t=0)
indices = np.argsort(magnitudes)[1:11]  # skip lowest (t=0)
pseudo_zeros = t_values[indices]

# Explicit formula correction (sum over pseudo-zeros)
def explicit_correction(n, pseudo_zeros, x_base=90):
    correction = 0.0
    for t in pseudo_zeros:
        rho = 0.5 + 1j * t
        x = x_base * n + 17  # for class 17
        li_rho = li_approx(x**rho)  # complex li, but approx via series
        correction -= (li_rho / rho).real  # take real part
    return correction

# Approximate pi_class(n) = main + explicit correction
pi_explicit = approx_main + np.array([explicit_correction(n, pseudo_zeros) for n in n_values])

# Plot or print
print("Sample: n=epoch-1, Explicit Approx:", pi_explicit[-1])

