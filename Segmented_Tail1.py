import math
import numpy as np
import matplotlib.pyplot as plt

K = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89]

def get_operators(k):
    # ... same as before ...
    # (I'll use the same 24 fixed operators for k=17 as in your code)

operators = [
    (72,-1,17),(72,-1,91),(108,29,19),(108,29,53),
    (72,11,37),(72,11,71),(18,0,73),(18,0,89),
    (102,20,11),(102,20,67),(138,52,13),(138,52,29),
    (102,28,31),(102,28,47),(48,3,49),(48,3,83),
    (78,8,23),(78,8,79),(132,45,7),(132,45,41),
    (78,16,43),(78,16,59),(42,4,61),(42,4,77)
]

def compute_tail_marking(h_max=200, track_redundancy=True):
    tail_density_history = []
    raw_density_history = []
    redundant_new_ops = []

    cumulative_lcm = 1                     # Running LCM of all previous periods
    previous_periods = []

    for h in range(1, h_max + 1):
        epoch = 90 * h * h - 12 * h + 1
        prev_epoch = 90 * (h-1)**2 - 12*(h-1) + 1 if h > 1 else 0
        tail_length = epoch - prev_epoch

        new_p = [z + 90*(h-1) for _,_,z in operators]   # the 24 new periods at this h

        # Count how many of the new operators are redundant (LCM-bound to previous)
        redundant = 0
        for p_new in new_p:
            if math.gcd(p_new, cumulative_lcm) > 1:
                redundant += 1
        redundant_new_ops.append(redundant)

        # Update cumulative LCM
        for p in new_p:
            cumulative_lcm = math.lcm(cumulative_lcm, p)
            previous_periods.append(p)

        # Effective marking density in tail ≈ sum 1/p over all previous + new non-redundant
        effective = sum(1/p for p in previous_periods) + sum(1/p for p in new_p if math.gcd(p, cumulative_lcm) == p)
        raw_density = min(effective, 10.0)   # cap for display

        tail_hole_prob = math.exp(-raw_density) if raw_density < 20 else 0

        tail_density_history.append(tail_hole_prob)
        raw_density_history.append(raw_density)

        print(f"h={h:3d} | Tail length={tail_length:6,d} | New redundant={redundant:2d}/24 | Raw density={raw_density:6.3f} | Hole prob={tail_hole_prob:.5f}")

    # Plot
    plt.figure(figsize=(12, 6))
    plt.plot(range(1, h_max+1), raw_density_history, 'b-', label='Raw Marking Density')
    plt.plot(range(1, h_max+1), tail_density_history, 'r-', label='Hole Probability in Tail')
    plt.axhline(4.36, color='b', linestyle='--', label='Asymptotic raw density ≈4.36')
    plt.xlabel('Epoch h'), plt.ylabel('Density / Probability')
    plt.title('Marking Pressure & Hole Probability in the Tail Epoch')
    plt.legend()
    plt.grid(True)
    plt.savefig('tail_marking_pressure.png')
    plt.show()

compute_tail_marking(h_max=150)