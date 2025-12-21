# Step 1: Import necessary libraries
# We need numpy for arrays and FFT, matplotlib for plotting (optional but useful for visualization)
import numpy as np
import matplotlib.pyplot as plt

# Step 2: Define the sieve parameters from the manuscripts
# From evolved_raindrops1.py and attachments, params for Class 17 and 19
# We'll use h=3 for a small epoch=775 to demonstrate; increase for larger datasets

h = 50  # Epoch parameter
epoch = 90 * (h * h) - 12 * h + 1  # 775 for h=3
limit = epoch

# Class 17 parameters (from evolved_raindrops1.py and TwinPrime_JULY.pdf Appendix)
params_17 = [
    (132, 45, 7), (102, 20, 11), (138, 52, 13), 
    (72, -1, 17), (108, 29, 19), (78, 8, 23), (138, 52, 29),
    (102, 28, 31), (72, 11, 37), (132, 45, 41),
    (78, 16, 43), (102, 28, 47),
    (48, 3, 49), (108, 29, 53), (78, 16, 59),
    (42, 4, 61), (102, 20, 67), (72, 11, 71),
    (18, 0, 73), (42, 4, 77), (78, 8, 79), (48, 3, 83),
    (18, 0, 89), (72, -1, 91) 
]

# Class 19 parameters
params_19 = [
    (70, -1, 91), (20, 0, 89), (74, 5, 83), (70, 7, 79), (56, 6, 77), (34, 3, 73), 
    (20, 0, 71), (106, 21, 67), (70, 13, 61), (110, 27, 59), (74, 15, 53), (70, 13, 49), 
    (56, 6, 47), (124, 40, 43), (110, 33, 41), (106, 31, 37), (70, 7, 31), (110, 33, 29), 
    (74, 5, 23), (70, -1, 19), (146, 59, 17), (124, 40, 13), (110, 27, 11), (106, 21, 7)
]

all_params = params_17 + params_19

# Step 3: Generate amplitude lists (markings) for Class 17 and 19
# This is the "dataset" - arrays where amplitude[n] = number of marks at address n
# Isolated lists per operator, then sum for each class

N_OPS = 48  # 24 per class
amplitude_lists = [[0] * limit for _ in range(N_OPS)]

def mark_chain(x, l, m, z, lst):
    y = 90 * (x * x) - l * x + m
    if y >= len(lst):
        return
    lst[int(y)] += 1
    p = z + 90 * (x - 1)
    n = 1
    while True:
        pos = y + p * n
        if pos >= len(lst):
            break
        lst[int(pos)] += 1
        n += 1

# Compute new_limit for x iterations (from quadratic root)
a, b, c = 90, -300, 250 - limit
new_limit = (-b + np.sqrt(b**2 - 4*a*c)) / (2*a)

# Run the sieve
for k in range(N_OPS):
    l, m, z = all_params[k]
    for x in range(1, int(new_limit.real) + 1):
        mark_chain(x, l, m, z, amplitude_lists[k])

# Sum for each class
class17_amp = [sum(amplitude_lists[k][n] for k in range(24)) for n in range(limit)]
class19_amp = [sum(amplitude_lists[k][n] for k in range(24, 48)) for n in range(limit)]

print(f"Generated amplitudes for epoch {epoch}:")
print(f"Class 17 mean amplitude: {np.mean(class17_amp):.4f}")
print(f"Class 19 mean amplitude: {np.mean(class19_amp):.4f}")

# Step 4: Perform Fourier Transform on the amplitude lists
# FFT decomposes the signal into frequencies; peaks show periodicities in markings
# Use numpy.fft for DFT

def compute_fft(signal):
    N = len(signal)
    Y = np.fft.fft(signal)  # DFT
    freq = np.fft.fftfreq(N, d=1)  # Frequencies (d=1 for address spacing)
    mag = np.abs(Y)  # Magnitudes
    # Sort non-DC frequencies by magnitude (exclude f=0)
    idx = np.argsort(mag)[::-1]  # Descending
    top_idx = [i for i in idx if freq[i] != 0][:10]  # Top 10 non-zero
    return freq[top_idx], mag[top_idx]

# Compute for Class 17
freq17, mag17 = compute_fft(class17_amp)
print("\nTop frequencies for Class 17 amplitudes:")
for f, m in zip(freq17, mag17):
    period = 1 / abs(f) if f != 0 else 'DC'
    print(f"Freq {f:.4f}: Mag {m:.2f} (period ≈{period:.2f})")

# Compute for Class 19
freq19, mag19 = compute_fft(class19_amp)
print("\nTop frequencies for Class 19 amplitudes:")
for f, m in zip(freq19, mag19):
    period = 1 / abs(f) if f != 0 else 'DC'
    print(f"Freq {f:.4f}: Mag {m:.2f} (period ≈{period:.2f})")

# Step 5: Optional - Plot the spectra for visualization
plt.figure(figsize=(12, 6))
plt.subplot(121)
plt.stem(freq17, mag17)
plt.title('Class 17 FFT Spectrum (Top 10)')
plt.xlabel('Frequency')
plt.ylabel('Magnitude')

plt.subplot(122)
plt.stem(freq19, mag19)
plt.title('Class 19 FFT Spectrum (Top 10)')
plt.xlabel('Frequency')
plt.ylabel('Magnitude')

plt.tight_layout()
plt.show()

# Explanation of steps:
# 1. Generate amplitudes: Run sieve to get marking counts per n (class17_amp, class19_amp).
# 2. FFT: np.fft.fft computes complex DFT; take abs for magnitude.
# 3. Frequencies: np.fft.fftfreq gives cycles per address; periods = 1/|freq|.
# 4. Top peaks: Sorted by mag, excluding DC (total sum), show dominant periodicities from operator chains.
# These align with 90/p fractions from base residues z.