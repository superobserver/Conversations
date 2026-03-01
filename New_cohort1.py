import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Operator 7 params (l,m from sample silo 17, z=7; generalize as needed)
l = 132  # From your ops for z=7
m = 45
z = 7

# Window params
n0 = 1e6  # Start
w = 1000  # Width
x_max = 1000  # Scan x

# Compute roll-in/out
insertion_rates = []
overlap_periods = []
remnants = []  # Prints in window
for x in range(1, x_max + 1):
    y = 90 * x * x - l * x + m
    p = z + 90 * (x - 1)
    if y > n0 + w:  # Roll-out (exclusion)
        continue
    # Insertion rate: fraction marked ~ w / p if intersects
    intersects = False
    prints = []
    for n in range(int(w / p) + 2):
        mark = y + n * p
        if n0 <= mark < n0 + w:
            intersects = True
            prints.append(mark - n0)
    if intersects:
        insertion_rates.append(len(prints) / w)
        # Overlap period with initial (x=1)
        p1 = z
        overlap = math.lcm(p, p1) if p and p1 else 0
        overlap_periods.append(overlap)
        remnants.append(prints)
    else:
        # Exclusion: distance from n-th primitive remnant
        p1 = z
        n_iter = math.floor((n0 - y) / p1)
        remnant = y + n_iter * p1
        dist = n0 - remnant
        print(f"x={x}: Excluded, Escape Dist ~{dist:.0f}")

print(f"Interacting x: {len(insertion_rates)} / {x_max-1}")
print("Sample Insertion Rates:", insertion_rates[:5])
print("Sample Overlap Periods:", overlap_periods[:5])
print("Sample Remnants (Prints):", remnants[0])

# Plot 1: Insertion Rate Decay
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(insertion_rates)+1), insertion_rates, 'b-o', label='Insertion Rate')
plt.xlabel('x'), plt.ylabel('Rate'), plt.title('Decay of Insertion Rate for Operator 7')
plt.legend(), plt.savefig('insertion_decay.png')

# Plot 2: Overlap Heatmap (Sample Remnants)
data = np.zeros((len(remnants), w))
for i, rem in enumerate(remnants):
    for pos in rem:
        data[i, pos] = 1
plt.figure(figsize=(12, 8))
sns.heatmap(data, cbar=False, cmap='binary')
plt.title('Remnant Prints in Window per Interacting x')
plt.xlabel('Window Index'), plt.ylabel('x Index')
plt.savefig('remnant_heatmap.png')

print("Plots Saved: 'insertion_decay.png' (rate decay), 'remnant_heatmap.png' (quadratic ties)") 