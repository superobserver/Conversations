import cmath
import math
import numpy as np
from scipy.optimize import curve_fit
import random

# Damped sine fit function
def damped_sine(x, a, b, c, d):
    return a * np.exp(-b * x) * np.sin(c * x + d)

# List of all 24 operator parameter tuples (l, m, z) from the sieve
operators = [
    (72, -1, 17), (72, -1, 91),      # 17,91
    (108, 29, 19), (108, 29, 53),    # 19,53
    (72, 11, 37), (72, 11, 71),      # 37,71
    (18, 0, 73), (18, 0, 89),        # 73,89
    (102, 20, 11), (102, 20, 67),    # 11,67
    (138, 52, 13), (138, 52, 29),    # 13,29
    (102, 28, 31), (102, 28, 47),    # 31,47
    (48, 3, 49), (48, 3, 83),        # 49,83
    (78, 8, 23), (78, 8, 79),        # 23,79
    (132, 45, 7), (132, 45, 41),     # 7,41
    (78, 16, 43), (78, 16, 59),      # 43,59
    (42, 4, 61), (42, 4, 77)         # 61,77
]

def compute_new_limit(h):
    epoch = 90 * (h * h) - 12 * h + 1
    limit = epoch
    a, b, c = 90, -300, 250 - limit
    d = (b**2) - (4 * a * c)
    sol2 = (-b + cmath.sqrt(d)) / (2 * a)
    return sol2.real, limit

def run_sieve(h):
    new_limit, limit = compute_new_limit(h)
    list17 = [0] * (int(limit) + 100)

    def drLD(x, l, m, z, listvar):
        y = 90 * (x * x) - l * x + m
        if 0 <= y < len(listvar):
            listvar[int(y)] += 1
        p = z + (90 * (x - 1))
        for n in range(1, int(((limit - y) / p) + 1) + 1):
            idx = y + (p * n)
            if 0 <= idx < len(listvar):
                listvar[int(idx)] += 1

    for x in range(1, int(new_limit) + 1):
        for l, m, z in operators:
            drLD(x, l, m, z, list17)

    list17 = list17[:int(limit)]
    holes = [i for i, amp in enumerate(list17) if amp == 0 and i > 1000]
    composites = [i for i, amp in enumerate(list17) if amp > 0 and i > 1000]
    return holes, composites, list17, new_limit

def analyze_oscillation(n, new_limit):
    op_results = []
    for op_idx, (l, m, z) in enumerate(operators):
        frac_devs = []
        for x in range(1, int(new_limit) + 1):
            y = 90 * (x * x) - l * x + m
            if y >= n:
                break
            p = z + 90 * (x - 1)
            if p == 0:  # Avoid div by zero, rare
                continue
            b = (n - y) / p
            if not b.is_integer():
                frac = b - math.floor(b)
                signed_dev = frac - 0.5
                frac_devs.append(signed_dev)

        if len(frac_devs) < 3:
            op_results.append((op_idx, 0, False, frac_devs))  # damp_rate, is_convergent, devs
            continue

        x_data = np.arange(len(frac_devs))
        try:
            popt, _ = curve_fit(damped_sine, x_data, np.array(frac_devs),
                                p0=[0.5, 0.01, 2 * np.pi / len(frac_devs), 0],
                                maxfev=10000, bounds=(-np.inf, np.inf))
            damp_rate = popt[1]  # Positive b indicates damping/convergence
            is_convergent = damp_rate > 0
        except:
            damp_rate = 0
            is_convergent = False

        op_results.append((op_idx, damp_rate, is_convergent, frac_devs))

    return op_results

def main():
    h = int(input("Enter h for sieve (e.g., 50 for epoch ~224k, 100 for ~900k): "))
    print(f"Running sieve for h={h}...")
    holes, composites, amplitudes, new_limit = run_sieve(h)
    print(f"Epoch: {90 * h**2 - 12 * h + 1}, Holes: {len(holes)}, Composites: {len(composites)}")

    # Sample 20 each for analysis
    sample_holes = random.sample(holes, min(20, len(holes)))
    sample_comps = random.sample(composites, min(20, len(composites)))

    print("\nHole Analysis:")
    for n in sample_holes:
        results = analyze_oscillation(n, new_limit)
        conv_count = sum(1 for _, _, conv, _ in results if conv)
        avg_damp = np.mean([d for _, d, _, _ in results if d != 0]) or 0
        print(f"n={n}: Conv ops={conv_count}/24, Avg damp={avg_damp:.4f}")

    print("\nComposite Analysis:")
    for n in sample_comps:
        results = analyze_oscillation(n, new_limit)
        conv_count = sum(1 for _, _, conv, _ in results if conv)
        avg_damp = np.mean([d for _, d, _, _ in results if d != 0]) or 0
        amp = amplitudes[n]  # Actual amplitude (sum of hits)
        print(f"n={n} (amp={amp}): Conv ops={conv_count}/24, Avg damp={avg_damp:.4f}")

if __name__ == "__main__":
    main()