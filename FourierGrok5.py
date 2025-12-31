import numpy as np
import matplotlib.pyplot as plt

# Series approximation for li(x) = ∫_2^x dt / ln(t) (offset version, accurate for x > 2)
def li_series(x, terms=20):
    if x <= 2:
        return 0.0
    logx = np.log(x)
    s = 0.0
    for k in range(terms):
        s += x / (logx ** (k + 1)) * np.math.factorial(k) / np.math.factorial(k + 1)
    return x / logx + s

# Vectorized version
li_approx = np.vectorize(li_series)

# Full sieve for class 17 (from your params_17 list)
def run_sieve(h=50):
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
        y = 90 * x * x - l * x + m
        if y >= len(lst): return
        lst[int(y)] += 1
        p = z + 90 * (x - 1)
        n = 1
        while y + p * n < len(lst):
            lst[int(y + p * n)] += 1
            n += 1
    
    for k, (l, m, z) in enumerate(params_17):
        for x in range(1, new_limit + 1):
            mark(x, l, m, z, amplitude_lists[k])
    
    class17_amp = np.array([sum(amplitude_lists[k][n] for k in range(24)) for n in range(epoch)])
    mean_amp = class17_amp.mean()
    print(f"Class 17 mean amplitude: {mean_amp:.4f}")
    
    return class17_amp, epoch

# Simplified continued_modular_zeta (from your FourierGrok3.py)
def continued_modular_zeta(signal, s_grid):
    N = len(signal)
    Y = np.fft.fft(signal - signal.mean())
    freq = np.fft.fftfreq(N)
    Z = np.zeros(len(s_grid), dtype=complex)
    for i, s in enumerate(s_grid):
        Z[i] = np.sum(Y * np.exp(-2j * np.pi * s * freq)) / N
    return Z + 1  # trivial term

# Main computation
h = 50  # or your larger h
class17_amp, epoch = run_sieve(h)

n_values = np.arange(epoch)
N_values = 90 * n_values + 17
approx_main = li_approx(N_values) / 24  # PNT in AP: density 1/φ(90)=1/24

# FFT for pseudo-zeros
s_grid = 0.5 + 1j * np.linspace(0, 120, 600)  # your grid
Z17 = continued_modular_zeta(class17_amp, s_grid)
magnitudes = np.abs(Z17)
t_values = s_grid.imag
indices = np.argsort(magnitudes)[1:11]  # top 10 minima, skip t=0
pseudo_zeros = t_values[indices]
print("Pseudo-zeros (imag parts):", pseudo_zeros)

# Explicit correction (simplified; real part only, no complex li yet)
def explicit_correction(n, pseudo_zeros):
    correction = 0.0
    for t in pseudo_zeros:
        # Approximate li(x^rho) ~ li(x^{0.5}) * cos(t ln x) / sqrt(x) + sin term, but use series
        x = 90 * n + 17
        rho_real = 0.5
        x_rho = x ** rho_real * np.exp(1j * t * np.log(x))  # complex
        li_rho = li_series(np.abs(x_rho))  # magnitude approx
        correction -= li_rho * np.cos(t * np.log(x)) / (rho_real + 1)  # crude real part
    return correction

pi_explicit = approx_main + np.array([explicit_correction(n, pseudo_zeros) for n in n_values])

# Exact count
primes17_exact = np.cumsum(class17_amp == 0)

# Results
print("Sample: n=1000, Exact:", primes17_exact[min(1000, epoch-1)],
      "Main:", approx_main[min(1000, epoch-1)],
      "Explicit Approx:", pi_explicit[min(1000, epoch-1)])
print("Sample: n=epoch-1, Exact:", primes17_exact[-1],
      "Main:", approx_main[-1],
      "Explicit Approx:", pi_explicit[-1])

# Plot (optional)
plt.figure(figsize=(12,6))
plt.plot(n_values, primes17_exact, label='Exact')
plt.plot(n_values, approx_main, '--', label='Main Term')
plt.plot(n_values, pi_explicit, '-.', label='Explicit Correction')
plt.legend()
plt.show()

import numpy as np

# Adjusted li(x) series (as before, vectorized)
def li_series(x, terms=20):
    if np.isscalar(x):
        x = np.array([x])
    result = np.zeros_like(x, dtype=float)
    mask = x > 2
    logx = np.log(x[mask])
    for k in range(terms):
        result[mask] += x[mask] / (logx ** (k + 1)) * np.math.factorial(k)
    return result / logx if np.isscalar(x) else result

# Run your sieve to get class17_amp, epoch (assume done)
# For illustration, use your values
epoch = 224401
primes17_exact = 53455  # from your output
mean_amp = 1.7048
hole_density = primes17_exact / epoch  # ~0.238

n_values = np.arange(epoch)
N_values = 90 * n_values + 17
approx_main_adjusted = li_series(N_values) * hole_density / 24  # scale by empirical density

# Pseudo-zeros from your output
pseudo_zeros = np.array([0.20033389, 0.40066778, 0.60100167, 0.80133556, 1.00166945,
                         1.20200334, 1.40233723, 1.60267112, 1.80300501, 2.0033389])

def explicit_correction_phase(n, pseudo_zeros):
    correction = 0.0
    x = 90 * n + 17
    if x <= 2: return 0.0
    sqrt_x = np.sqrt(x)
    ln_x = np.log(x)
    for t in pseudo_zeros:
        phase = t * ln_x
        amp = sqrt_x * np.cos(phase) / ln_x
        rho = 0.5 + 1j * t
        correction -= amp / abs(rho)  # leading oscillatory term
    return correction

pi_explicit_adjusted = approx_main_adjusted + np.array([explicit_correction_phase(n, pseudo_zeros) for n in n_values])

# Results (run locally for full array)
print("Adjusted Main at epoch-1:", approx_main_adjusted[-1])
print("Explicit Approx at epoch-1:", pi_explicit_adjusted[-1])
print("Exact:", primes17_exact)