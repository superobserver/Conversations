import math
import numpy as np

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

# Simulate race for max_x=10 (small for speed; adjust as needed)
max_x = 90
hole_history = {k: [0] for k in ALL_CLASSES}  # Start with 0 holes at x=0
for h in range(1, max_x + 1):
    for k in ALL_CLASSES:
        cumulative_holes = compute_holes(h, k)
        hole_history[k].append(cumulative_holes)

# Derive positions
position_counts = {k: {pos: 0 for pos in range(1, 25)} for k in ALL_CLASSES}
for x in range(1, max_x + 1):
    counts = {k: hole_history[k][x] for k in ALL_CLASSES}
    sorted_k = sorted(counts, key=counts.get, reverse=True)
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

# Compute leaderboard: For each position, silo with most time there
leaderboard = {}
for pos in range(1, 25):
    max_time = 0
    leader_silo = None
    for k in ALL_CLASSES:
        time = position_counts[k][pos]
        if time > max_time:
            max_time = time
            leader_silo = k
    leaderboard[pos] = (leader_silo, max_time)

# Output leaderboard
print("Leaderboard:")
for pos, (silo, time) in leaderboard.items():
    print(f"Position {pos}: Silo {silo} ({time} iterations)")

# Conjecture analysis (simple check for preferences)
print("\nConjecture Check: Spectral Preferences")
for k in ALL_CLASSES:
    prefs = sorted(position_counts[k].items(), key=lambda item: item[1], reverse=True)[:3]  # Top 3 positions
    print(f"Silo {k}: Top Positions {prefs}")