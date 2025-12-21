import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from sympy import isprime

# Parameters
h = 10  # Epoch scale, adjust for larger range
epoch = 90 * h**2 - 12 * h + 1
n_range = np.arange(101)  # Simulate up to n=100
residues = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]
k = 11  # Class for A201804
k_twin = 13  # Twin class for A224854

# Quadratic operator pairs for k=11 (example pairs from CRT-derived r,s)
pairs_k = [(1, 11), (7, 41), (13, 29)]  # Subset for illustration, expand as needed
pairs_k_twin = [(1, 13), (7, 43), (11, 37)]  # For k=13

# Sine wave model: ψ(n) = A sin(2π n / p + φ), rectified to +1/0
def sine_operator(n, p, y_start, A=1):
    phase = 2 * np.pi * (n - y_start) / p
    return np.where(np.sin(phase) > 0, A, 0)  # Rectify positive peaks to mark

# Amplitude sum over operators
def compute_amplitude(n, pairs, h):
    amp = np.zeros_like(n, dtype=float)
    x_max = int(h + 5/3)
    for x in range(1, x_max + 1):
        for r, s in pairs:
            l = 180 - (r + s)
            m = (r * s - k) // 90 + 90 - l  # CRT-adjusted start
            y = 90 * x**2 - l * x + m
            p = r + 90 * (x - 1)  # Approximate period from operator
            if y >= 0 and y < len(n):
                amp = np.maximum(amp, sine_operator(n, p, y % p))
    return amp

# Initialize animation
fig, ax = plt.subplots(figsize=(10, 6))
line, = ax.plot([], [], 'b-', label='Amplitude Sum')
ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
ax.set_xlim(0, 100)
ax.set_ylim(-0.1, 2.1)
ax.set_xlabel('Index n')
ax.set_ylabel('Amplitude')
ax.set_title('Sine Operator Sum for Class 11 and Twin 13')
ax.legend()

# Animation update function
def update(frame):
    if frame < len(pairs_k):
        amp_k = compute_amplitude(n_range, pairs_k[:frame + 1], h)
        amp_twin = compute_amplitude(n_range, pairs_k_twin[:frame + 1], h)
        total_amp = np.maximum(amp_k, amp_twin)  # Union of markings
    else:
        amp_k = compute_amplitude(n_range, pairs_k, h)
        amp_twin = compute_amplitude(n_range, pairs_k_twin, h)
        total_amp = np.maximum(amp_k, amp_twin)
    line.set_data(n_range, total_amp)
    return line,

# Animate the buildup
ani = FuncAnimation(fig, update, frames=len(pairs_k) + 1, interval=500, blit=True)
ani.save('sine_operator_sums.gif', writer='pillow')
print("Animation saved as sine_operator_sums.gif")

# Static plot with hole verification
plt.figure(figsize=(10, 6))
amp_k = compute_amplitude(n_range, pairs_k, h)
amp_twin = compute_amplitude(n_range, pairs_k_twin, h)
total_amp = np.maximum(amp_k, amp_twin)
plt.plot(n_range, total_amp, 'b-', label='Amplitude Sum')
for n in n_range:
    if isprime(90 * n + 11) and isprime(90 * n + 13):
        plt.plot(n, total_amp[n], 'ro', label='Twin Hole' if n == 0 else "")
plt.axhline(y=0, color='k', linestyle='-', alpha=0.3)
plt.xlabel('Index n')
plt.ylabel('Amplitude')
plt.title('Sine Operator Sum for Classes 11 and 13 with Twin Holes')
plt.legend()
plt.savefig('sine_operator_sums.png')
plt.show()
print("Static plot saved as sine_operator_sums.png")

# 3D surface for geometric insight
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
X, Y = np.meshgrid(n_range, np.arange(len(pairs_k) + len(pairs_k_twin)))
Z = np.zeros_like(X)
for i, (r, s) in enumerate(pairs_k + pairs_k_twin):
    l = 180 - (r + s)
    m = (r * s - k) // 90 + 90 - l
    for x in range(1, h + 1):
        y = 90 * x**2 - l * x + m
        if 0 <= y < len(n_range):
            Z[i, y] = 1  # Mark operator contribution
ax.plot_surface(X, Y, Z, cmap='viridis')
ax.set_xlabel('Index n')
ax.set_ylabel('Operator Pair Index')
ax.set_zlabel('Amplitude Contribution')
ax.set_title('3D Geometry of Sine Operator Pairs')
plt.show()