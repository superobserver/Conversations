import math
import numpy as np

# Define the 24 classes
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

def compute_cumulative_holes_up_to_h(max_h):
    # For each h from 1 to max_h, compute hole counts for each class
    # But to get per "x+1" (per operator level), we need to think in terms of x iterations
    # Since epochs grow with h, but to simulate "race" per x+1, let's iterate x from 1 to some max_x
    # For each x, compute incremental holes added by that x's operators across classes
    # But since holes are unmarked, it's cumulative unmarked up to current epoch approx

    # Choose max_x such that computation is feasible (e.g., 50, as before)
    max_x = 50
    # For each x, the "epoch" effectively covered is up to approx 90 x^2
    # But to track cumulative holes, we need to run the sieve up to each x's contribution

    # To simulate the "race", we'll run the sieve for increasing max_x, computing total holes for each class at each step
    # This is compute-intensive, but for max_x=50, epoch~112k, 24 classes, feasible

    # Dict to store cumulative holes per class per max_x
    cumulative_holes = {k: [] for k in ALL_CLASSES}
    spreads = []  # max - min holes per max_x
    max_gaps = []  # largest gap between consecutive ranked positions
    rel_spreads = []  # spread / avg_holes

    for current_x in range(1, max_x + 1):
        hole_counts = []
        for k in ALL_CLASSES:
            # Compute holes with new_limit = current_x
            epoch = 90 * current_x * current_x  # Approx, but use full for accuracy
            amplitude = np.zeros(int(epoch + 100), dtype=int)
            ops = get_operators(k)
            for l, m, z in ops:
                for x in range(1, current_x + 1):
                    y = 90 * x * x - l * x + m
                    if 0 <= y < len(amplitude):
                        amplitude[int(y)] += 1
                    p = z + 90 * (x - 1)
                    if p <= 0: continue
                    n = 1
                    yy = y + p * n
                    while yy < len(amplitude):
                        amplitude[int(yy)] += 1
                        n += 1
                        yy = y + p * n
            holes = np.sum(amplitude == 0)
            hole_counts.append(holes)
        cumulative_holes[current_x] = hole_counts  # Store per x
        
        # Compute rankings and metrics
        sorted_holes = sorted(hole_counts, reverse=True)
        spread = sorted_holes[0] - sorted_holes[-1]
        spreads.append(spread)
        
        gaps = [sorted_holes[i] - sorted_holes[i+1] for i in range(len(sorted_holes)-1)]
        max_gap = max(gaps) if gaps else 0
        max_gaps.append(max_gap)
        
        avg_holes = np.mean(hole_counts)
        rel_spread = spread / avg_holes if avg_holes > 0 else 0
        rel_spreads.append(rel_spread)

    # Output results
    print("Spreads (max - min):", spreads)
    print("Max Gaps:", max_gaps)
    print("Relative Spreads:", rel_spreads)
    
    # Plot
    xs = list(range(1, max_x + 1))
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(xs, spreads, label='Spread (max-min)')
    plt.plot(xs, max_gaps, label='Max Gap')
    plt.plot(xs, rel_spreads, label='Rel Spread')
    plt.xlabel('x iteration')
    plt.ylabel('Value')
    plt.legend()
    plt.savefig('race_metrics.png')
    plt.close()
    
    # For positions: Track rank changes, but for brevity, print for last x
    last_holes = hole_counts
    ranked = sorted(zip(ALL_CLASSES, last_holes), key=lambda p: p[1], reverse=True)
    print("Final Rankings:", ranked)
    
    return "Done"
    
    
compute_cumulative_holes_up_to_h(50)