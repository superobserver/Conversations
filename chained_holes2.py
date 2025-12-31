import numpy as np
import matplotlib.pyplot as plt

# Operator params for class 17 (from code snippets)
params_17 = [
    (132,45,7),(102,20,11),(138,52,13),(72,-1,17),(108,29,19),(78,8,23),
    (138,52,29),(102,28,31),(72,11,37),(132,45,41),(78,16,43),(102,28,47),
    (48,3,49),(108,29,53),(78,16,59),(42,4,61),(102,20,67),(72,11,71),
    (18,0,73),(42,4,77),(78,8,79),(48,3,83),(18,0,89),(72,-1,91)
]

# Simulate small epoch for visualization (E=10000, n=0 to 9999)
E = 10000
amp = np.zeros(E)

def mark_sim(n, l, m, z):
    # Simplified mark: add mark if n in chain (y + k*p <= E, p~z)
    p = z  # base period approximation
    y = 90 * (n // 10)**2 - l * (n // 10) + m  # mock quadratic
    if y >= E: return
    pos = int(y)
    while pos < E:
        amp[pos] += 1
        pos += p

for l, m, z in params_17:
    for n in range(E):
        mark_sim(n, l, m, z)

mean_amp = amp.mean()

# FFT of de-meaned signal
signal_demean = amp - mean_amp
Y = np.fft.fft(signal_demean)
freq = np.fft.fftfreq(E)

# Magnitude spectrum (positive frequencies only)
mag = np.abs(Y[:E//2])
f_pos = freq[:E//2]

# Chain peaks (high mag > 5*mean)
peak_threshold = 5 * np.mean(mag)
peak_freq = f_pos[mag > peak_threshold]

# Anti-spectrum (low mag < mean/5)
anti_threshold = np.mean(mag) / 5
anti_freq = f_pos[mag < anti_threshold]

# Prime-rich indices (n where signal low, i.e., anti-spectrum alignment)
# Approximate: indices with low local amplitude
low_amp_mask = amp < mean_amp / 2  # heuristic for prime-rich
prime_rich_n = np.where(low_amp_mask)[0]

print(f"Number of peak frequencies (chain peaks): {len(peak_freq)}")
print(f"Number of anti-spectrum frequencies: {len(anti_freq)}")
print(f"Approximate prime-rich indices: {len(prime_rich_n)}")

# Visualization
fig, axs = plt.subplots(3, 1, figsize=(12, 10))

# 1. Amplitude signal
axs[0].plot(amp, label='Amplitude (Marks)', color='blue', alpha=0.7)
axs[0].axhline(mean_amp, color='red', linestyle='--', label='Mean')
axs[0].set_title('Amplitude Signal (Composite Marks)')
axs[0].set_xlabel('Index n')
axs[0].set_ylabel('Marks')
axs[0].legend()

# 2. FFT Spectrum
axs[1].plot(f_pos, mag, label='Magnitude Spectrum', color='green')
axs[1].scatter(peak_freq, mag[mag > peak_threshold], color='orange', label='Chain Peaks')
axs[1].scatter(anti_freq, mag[mag < anti_threshold], color='purple', alpha=0.5, label='Anti-Spectrum')
axs[1].set_title('Fourier Spectrum: Chain Peaks vs Anti-Spectrum')
axs[1].set_xlabel('Frequency f')
axs[1].set_ylabel('|Y(f)|')
axs[1].legend()

# 3. Prime-rich indices
axs[2].scatter(prime_rich_n, amp[prime_rich_n], color='red', s=10, label='Prime-Rich n (Low Amp)')
axs[2].plot(amp, color='gray', alpha=0.3)
axs[2].set_title('Prime-Rich Indices (Anti-Spectrum Alignment)')
axs[2].set_xlabel('Index n')
axs[2].set_ylabel('Amplitude')
axs[2].legend()

plt.tight_layout()
plt.show()