import math
import numpy as np
import matplotlib.pyplot as plt

# Fixed operator pool
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

# Compute y for operator z=7 across siloes
def compute_y_for_op7(x_max=50):
    y_data = {}
    for k in ALL_CLASSES:
        ops = get_operators(k)
        # Find operators where z==7 or z==91 (if adjusted), but since z_eff for 1 is 91, for 7 it's 7
        relevant_ops = [op for op in ops if op[2] == 7]
        if not relevant_ops:
            y_data[k] = []
            continue
        # Assume one per silo for simplicity; take first if multiple
        l, m, _ = relevant_ops[0]
        ys = []
        for x in range(1, x_max + 1):
            y = 90 * x * x - l * x + m
            ys.append(y)
        y_data[k] = ys
    return y_data

# Run computation
y_data = compute_y_for_op7(x_max=50)

# Print summary
print("Distribution of y for operator 7 across siloes:")
for k, ys in y_data.items():
    if ys:
        min_y = min(ys)
        max_y = max(ys)
        avg_y = np.mean(ys)
        print(f"Silo {k}: Min y={min_y}, Max y={max_y}, Avg y={avg_y:.2f}")
    else:
        print(f"Silo {k}: Operator 7 not present")

# Plot 24 line graphs (or fewer if not present in all)
fig, axs = plt.subplots(6, 4, figsize=(20, 30))  # 6x4 grid for 24 plots
axs = axs.flatten()
x_vals = list(range(1, 51))  # x=1 to 50
plot_idx = 0
for k, ys in y_data.items():
    if ys:
        axs[plot_idx].plot(x_vals, ys, marker='o', linestyle='-', label=f'Silo {k}')
        axs[plot_idx].set_title(f'y vs x for Op 7 in Silo {k}')
        axs[plot_idx].set_xlabel('x')
        axs[plot_idx].set_ylabel('y')
        axs[plot_idx].legend()
        axs[plot_idx].grid(True)
        plot_idx += 1

plt.tight_layout()
plt.savefig('op7_y_distributions.png')
print("\nGraph saved as 'op7_y_distributions.png'. Description: 24 line plots (or fewer if op7 absent), each showing y(x) = 90x² - l x + m for the (l,m) associated with z=7 in that silo. Nodes at integer x+1 (but plotted at x for simplicity). Lines form quadratic curves, bounded in a 'cone' shape diverging from origin.")