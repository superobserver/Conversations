import math
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Fixed pool (24 classes coprime to 90)
ALL_CLASSES = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]

# Mock incremental hole computation (replace with your real sieve function)
def compute_holes_incremental(x, k, prev_holes):
    base = int(3.75 * x**2 / (9 + 2 * math.log(x+1)) / 24)
    fluctuation = int(20 * math.sin(0.3 * x + k)) + np.random.randint(-15, 15)
    return prev_holes + max(0, base + fluctuation)

# Simulate race
def run_race(max_x=200):
    hole_history = {k: [0] for k in ALL_CLASSES}  # holes at x=0
    position_history = {k: [] for k in ALL_CLASSES}  # rank at each x

    for x in range(1, max_x + 1):
        counts = {}
        for k in ALL_CLASSES:
            new_holes = compute_holes_incremental(x, k, hole_history[k][-1])
            hole_history[k].append(new_holes)
            counts[k] = new_holes

        sorted_k = sorted(counts, key=counts.get, reverse=True)
        rank_dict = {}
        current_rank = 1
        prev_count = None
        for idx, k in enumerate(sorted_k):
            count = counts[k]
            if count != prev_count:
                current_rank = idx + 1
            rank_dict[k] = current_rank
            prev_count = count

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

    pos_leader = {}
    for pos in range(1, 25):
        if pos_leaderboard[pos]:
            best_k = max(pos_leaderboard[pos], key=pos_leaderboard[pos].get)
            pos_leader[pos] = (best_k, pos_leaderboard[pos][best_k])

    return pos_leader, silo_preferences, pos_leaderboard

# Visualization: 24 Bar Graphs + Summary Heatmap
def plot_24_bars(pos_leaderboard, max_x):
    for k in ALL_CLASSES:
        fig, ax = plt.subplots(figsize=(10, 5))
        counts = [pos_leaderboard[pos].get(k, 0) for pos in range(1, 25)]
        ax.bar(range(1, 25), counts, color='steelblue')
        ax.set_title(f'Silo {k}: Time Spent in Each Position (over {max_x} iterations)')
        ax.set_xlabel('Position (1 = Leader, 24 = Last)')
        ax.set_ylabel('Number of Iterations')
        ax.set_xticks(range(1, 25))
        plt.tight_layout()
        plt.savefig(f'silo_{k}_position_bar.png')
        plt.close(fig)
    print(f"Generated 24 bar graphs: 'silo_k_position_bar.png' for each k in {ALL_CLASSES}")

    # Summary Heatmap
    data = np.array([[pos_leaderboard[pos].get(k, 0) for pos in range(1, 25)] for k in ALL_CLASSES])
    plt.figure(figsize=(14, 10))
    sns.heatmap(data, xticklabels=range(1,25), yticklabels=ALL_CLASSES, cmap='YlGnBu', annot=False)
    plt.title("Heatmap: Time Spent in Positions by Silo")
    plt.xlabel("Position (1=Leader, 24=Last)")
    plt.ylabel("Silo")
    plt.savefig('position_heatmap_summary.png')
    plt.close()

# Run everything
print("Running Race Simulation (x=1 to 200)...")
position_history = run_race(max_x=200)

print("\nLeaderboard: Position → (Silo, Time Spent)")
pos_leader, silo_pref, pos_leaderboard = analyze_preferences(position_history)
for pos in range(1, 25):
    silo, time = pos_leader[pos]
    print(f"Position {pos}: Silo {silo} ({time} iterations)")

print("\nConjecture Check: Spectral Preferences per Silo (Top 3 Positions)")
for k in ALL_CLASSES:
    print(f"Silo {k}: Top Positions {silo_pref[k]}")

plot_24_bars(pos_leaderboard, 200)
print("\nPlots Saved:")
print("- 24 individual bar graphs: 'silo_k_position_bar.png' for each silo")
print("- 'position_heatmap_summary.png': Overview of all siloes and positions")

# Sample LCM Matrix for Operators at x=10 (Silos 11 and 13)
K = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89]
p_11 = [z + 90*9 for z in K]  # p(x=10) for silo 11 (approx)
p_13 = [z + 90*9 for z in K]  # Similar, but variances in effective marking
lcm_matrix = [[math.lcm(p1, p2) for p2 in p_13] for p1 in p_11]
print("Sample LCM Matrix Variations (subset):", lcm_matrix[0][:5])  # Off-diagonals show density mods