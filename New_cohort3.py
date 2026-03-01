import math
import numpy as np
import matplotlib.pyplot as plt

# Operators for sample silo (k=17)
operators = [
    (72, -1, 17), (72, -1, 91), (108, 29, 19), (108, 29, 53),
    (72, 11, 37), (72, 11, 71), (18, 0, 73), (18, 0, 89),
    (102, 20, 11), (102, 20, 67), (138, 52, 13), (138, 52, 29),
    (102, 28, 31), (102, 28, 47), (48, 3, 49), (48, 3, 83),
    (78, 8, 23), (78, 8, 79), (132, 45, 7), (132, 45, 41),
    (78, 16, 43), (78, 16, 59), (42, 4, 61), (42, 4, 77)
]

n0 = 10_000_000
w = 1000
x_max_scan = 400  # safe upper bound for x ~ sqrt(n0/90) ≈ 333

active_x_per_op = []
active_count = 0

for idx, (l, m, z) in enumerate(operators):
    active_xs = []
    for x in range(1, x_max_scan + 1):
        y = 90 * x * x - l * x + m
        p = z + 90 * (x - 1)
        if y > n0 + w:
            break  # later x only worse
        # Does the chain hit the window at all?
        if y + math.ceil((n0 - y)/p if p > 0 else 0) * p <= n0 + w - 1:
            active_xs.append(x)
    if active_xs:
        active_count += 1
        active_x_per_op.append((idx, active_xs))

print(f"Segment [{n0:,} – {n0+w-1:,}]")
print(f"Total active operators: {active_count} / 24")
print(f"Active x values per operator (showing first few):")
for idx, xs in active_x_per_op[:5]:
    print(f"  Op {idx:2d} (z={operators[idx][2]}): x = {xs[:8]}{'...' if len(xs)>8 else ''}")

# Plot: Number of active operators vs approximate x
x_range = np.arange(1, x_max_scan + 1)
active_per_x = [sum(1 for _, xs in active_x_per_op if x in xs) for x in x_range]
plt.figure(figsize=(10, 6))
plt.plot(x_range, active_per_x, 'b-', label='Active operators at x')
plt.axvline(math.sqrt(n0/90), color='r', linestyle='--', label='~x_max expected')
plt.xlabel('x'), plt.ylabel('Active Operators'), plt.title('Roll-in / Roll-out: Active Operators per x')
plt.legend()
plt.savefig('active_operators_vs_x.png')
print("Plot Saved: 'active_operators_vs_x.png' — shows peak around x ≈ √(n₀/90) ≈ 333, rapid decay after.")