import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Fixed pool (24 classes coprime to 90)
ALL_CLASSES = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]

# Placeholder for hole computation (use your full sieve function here)
# For demonstration we use a mock that gives realistic variation
def compute_holes_incremental(x, k, prev_holes):
    # Mock: base PNT + modular fluctuation + noise
    base = int(3.75 * x**2 / (9 + 2 * math.log(x+1)) / 24)
    fluctuation = int(20 * math.sin(0.3 * x + k)) + np.random.randint(-15, 15)
    return prev_holes + max(0, base + fluctuation)

# Simulate race
def run_race(max_x=400):
    hole_history = {k: [0] for k in ALL_CLASSES}  # holes at x=0
    position_history = {k: [] for k in ALL_CLASSES}  # rank at each x

    for x in range(1, max_x + 1):
        counts = {}
        for k in ALL_CLASSES:
            new_holes = compute_holes_incremental(x, k, hole_history[k][-1])
            hole_history[k].append(new_holes)
            counts[k] = new_holes

        # Rank: 1=leader (most holes), 24=last
        sorted_k = sorted(counts, key=counts.get, reverse=True)
        rank_dict = {k: idx + 1 for idx, k in enumerate(sorted_k)}
        for k in ALL_CLASSES:
            position_history[k].append(rank_dict[k])

    return position_history

# Analyze leaderboard and preferences
def analyze_preferences(position_history):
    pos_leaderboard = {pos: {} for pos in range(1, 25)}
    silo_preferences = {}
    
    for k in ALL_CLASSES:
        counts = np.bincount(position_history[k], minlength=25)[1:]
        top3_idx = np.argsort(counts)[-3:][::-1] + 1
        top3_vals = counts[top3_idx - 1]
        silo_preferences[k] = list(zip(top3_idx, top3_vals))
        
        for pos, cnt in enumerate(counts, 1):
            pos_leaderboard[pos][k] = cnt

    # Position leaderboard: silo with most time in each position
    pos_leader = {}
    for pos in range(1, 25):
        if pos_leaderboard[pos]:
            best_k = max(pos_leaderboard[pos], key=pos_leaderboard[pos].get)
            pos_leader[pos] = (best_k, pos_leaderboard[pos][best_k])

    return pos_leader, silo_preferences

# Visualization
def plot_preferences(pos_leader, silo_preferences, max_x):
    # 1. Leaderboard heatmap
    plt.figure(figsize=(12, 8))
    data = np.zeros((24, 24))
    for pos in range(1, 25):
        for k_idx, k in enumerate(ALL_CLASSES):
            data[pos-1, k_idx] = pos_leader.get(pos, (0,0))[1] if pos_leader.get(pos, (0,0))[0] == k else 0
    sns.heatmap(data, xticklabels=ALL_CLASSES, yticklabels=range(1,25), cmap='YlOrRd')
    plt.title("Time Spent in Each Position by Silo (Leaderboard View)")
    plt.xlabel("Silo"), plt.ylabel("Position (1=Leader)"), plt.savefig('leaderboard_heatmap.png')
    plt.close()

    # 2. Positional entropy (lower = more spectral preference)
    entropies = {}
    for k in ALL_CLASSES:
        hist = np.bincount(position_history[k], minlength=25)[1:]
        hist = hist / hist.sum()
        ent = -np.sum(hist * np.log(hist + 1e-10))
        entropies[k] = ent

    plt.figure(figsize=(10, 6))
    plt.bar(ALL_CLASSES, [entropies[k] for k in ALL_CLASSES])
    plt.axhline(math.log(24), color='r', linestyle='--', label='Uniform Entropy (ln 24 ≈ 3.18)')
    plt.xlabel("Silo"), plt.ylabel("Positional Entropy"), plt.title("Spectral Preference: Lower Entropy = Stronger Preference")
    plt.legend(), plt.savefig('entropy_bar.png')
    plt.close()

# Run
print("Running the Silo Race Simulation...")
position_history = run_race(max_x=200)

print("\nLeaderboard: Position → (Silo, Time Spent)")
pos_leader, silo_pref = analyze_preferences(position_history)
for pos in range(1, 25):
    silo, time = pos_leader[pos]
    print(f"Position {pos}: Silo {silo} ({time} iterations)")

print("\nConjecture Check: Spectral Preferences per Silo")
for k in ALL_CLASSES:
    print(f"Silo {k}: Top Positions {silo_pref[k]}")

plot_preferences(pos_leader, silo_pref, 100)
print("\nPlots saved:")
print("- leaderboard_heatmap.png : Shows which silo dominates each position")
print("- entropy_bar.png : Lower bars indicate stronger spectral preference")