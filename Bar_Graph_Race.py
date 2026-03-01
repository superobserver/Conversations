import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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

# Simulate for x_max=50
x_max = 50
hole_history = {k: [0] for k in ALL_CLASSES}  # Start with 0 holes at x=0
for h in range(1, x_max + 1):
    for k in ALL_CLASSES:
        new_holes = compute_holes(h, k) - hole_history[k][-1]  # Incremental
        hole_history[k].append(hole_history[k][-1] + new_holes)

# Animation Setup
fig, ax = plt.subplots(figsize=(12, 8))
def animate(frame):
    ax.clear()
    counts = [hole_history[k][frame] for k in ALL_CLASSES]
    sorted_idx = np.argsort(counts)[::-1]
    sorted_classes = [ALL_CLASSES[i] for i in sorted_idx]
    sorted_counts = [counts[i] for i in sorted_idx]
    ax.bar(range(24), sorted_counts, tick_label=sorted_classes)
    ax.set_title(f'x={frame+1}: Hole Count Race (Reshuffled Bars)')
    ax.set_xlabel('Ranked Siloes'), ax.set_ylabel('Hole Count')
    ax.tick_params(axis='x', rotation=45)
    spread = max(sorted_counts) - min(sorted_counts)
    max_gap = max([sorted_counts[i] - sorted_counts[i+1] for i in range(23)])
    print(f"Frame x={frame+1}: Max Spread={spread}, Largest Gap={max_gap}")

ani = animation.FuncAnimation(fig, animate, frames=x_max, interval=200, repeat=False)
ani.save('race_reshuffle.gif', writer='pillow')
plt.close(fig)

print("Animation Saved: 'race_reshuffle.gif'—50 frames of reshuffling bars, labels by class.") 
print("Terminal Log: Per-frame metrics printed above—shows patterned fluctuations from L/24 influx.")