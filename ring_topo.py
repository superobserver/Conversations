import math
import numpy as np
import matplotlib.pyplot as plt

# Fixed pool
ALL_CLASSES = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]

def get_operators(k):
    operators = []
    seen = set()
    for z in ALL_CLASSES:
        try:
            o = (k * pow(z, -1, 90)) % 90
            if o not in ALL_CLASSES: continue
            pair = tuple(sorted([z, o]))
            if pair in seen: continue
            seen.add(pair)
            z_eff = 91 if z == 1 else z
            o_eff = 91 if o == 1 else o
            l = 180 - (z_eff + o_eff)
            m = 90 - (z_eff + o_eff) + (z_eff * o_eff - k) // 90
            operators.append((l, m, z_eff))
            if z != o:
                operators.append((l, m, o_eff))
        except ValueError:
            continue
    return operators

def compute_holes(h, k):
    epoch = 90 * h * h - 12 * h + 1
    a = 90
    b = -300
    c = 250 - epoch
    d = (b**2) - (4 * a * c)
    sol2 = (-b + math.sqrt(d)) / (2 * a) if d >= 0 else 0
    new_limit = sol2

    amplitude = [0] * int(epoch + 100)

    def drLD(x, l, m, z, listvar):
        y = 90 * x * x - l * x + m
        if 0 <= y < len(listvar):
            listvar[int(y)] += 1
        p = z + 90 * (x - 1)
        for n in range(1, int((epoch - y) / p) + 1):
            yy = y + p * n
            if 0 <= yy < len(listvar):
                listvar[int(yy)] += 1

    ops = get_operators(k)
    for l, m, z in ops:
        for x in range(1, int(new_limit) + 1):
            drLD(x, l, m, z, amplitude)

    amplitude = amplitude[:int(epoch)]
    holes = sum(1 for amp in amplitude if amp == 0)
    return holes

# Simulate race for max_x=10
max_x = 10
hole_history = {k: [] for k in ALL_CLASSES}
for h in range(1, max_x + 1):
    for k in ALL_CLASSES:
        holes = compute_holes(h, k)
        hole_history[k].append(holes)

# Derive time series and positions
first_place = []
last_place = []
position_counts = {k: {pos: 0 for pos in range(1, 25)} for k in ALL_CLASSES}
for x in range(max_x):
    counts = {k: hole_history[k][x] for k in ALL_CLASSES}
    sorted_k = sorted(counts, key=counts.get, reverse=True)
    first_place.append(sorted_k[0])
    last_place.append(sorted_k[-1])
    
    rank_dict = {}
    current_rank = 1
    prev_count = None
    for idx, k in enumerate(sorted_k):
        count = counts[k]
        if count != prev_count:
            current_rank = idx + 1
        rank_dict[k] = current_rank
        position_counts[k][current_rank] += 1
        prev_count = count

# Output
print("First Place Oscillation:", first_place)
print("Last Place Oscillation:", last_place)
print("Position Distribution for Class 7:", position_counts[7])

# Plot 1: Oscillation Time Series
plt.figure(figsize=(10, 6))
plt.plot(range(1, max_x+1), first_place, 'b-o', label='First Place Silo')
plt.plot(range(1, max_x+1), last_place, 'r-o', label='Last Place Silo')
plt.xlabel('x+1 Iteration'), plt.ylabel('Silo ID'), plt.legend(), plt.title('Time Series Oscillation: First/Last Positions')
plt.savefig('oscillation_series.png')

# Plot 2: Position Heatmap for All Siloes
positions = np.array([[position_counts[k][pos] for pos in range(1, 25)] for k in ALL_CLASSES])
plt.figure(figsize=(12, 8))
plt.imshow(positions, cmap='hot')
plt.colorbar(label='Time Spent')
plt.xticks(range(24), range(1, 25))
plt.yticks(range(24), ALL_CLASSES)
plt.title('Position Distribution Heatmap Across Siloes')
plt.savefig('position_heatmap.png')

print("Plots Saved: 'oscillation_series.png' and 'position_heatmap.png'")