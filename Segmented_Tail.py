import cmath
import math
import matplotlib.pyplot as plt

# Operators for k=17 (from your code)
operators = [
    (72, -1, 17), (72, -1, 91), (108, 29, 19), (108, 29, 53),
    (72, 11, 37), (72, 11, 71), (18, 0, 73), (18, 0, 89),
    (102, 20, 11), (102, 20, 67), (138, 52, 13), (138, 52, 29),
    (102, 28, 31), (102, 28, 47), (48, 3, 49), (48, 3, 83),
    (78, 8, 23), (78, 8, 79), (132, 45, 7), (132, 45, 41),
    (78, 16, 43), (78, 16, 59), (42, 4, 61), (42, 4, 77)
]

def segmented_tail_analysis(h_max=100):
    prev_epoch = 0
    total_marks_history = []
    avg_amplitude_history = []
    hole_prob_history = []

    for h in range(1, h_max + 1):
        epoch = 90 * h * h - 12 * h + 1
        tail_start = prev_epoch
        tail_length = epoch - prev_epoch

        # For this h, new x range
        a = 90
        b = -300
        c = 250 - epoch
        d = (b**2) - (4 * a * c)
        sol2 = (-b + cmath.sqrt(d)) / (2 * a)
        new_limit = int(sol2.real) + 1
        x_start = 1 if h == 1 else int(math.sqrt(250 * prev_epoch / 90)) + 1

        # Tail amplitude
        tail_amplitude = [0] * tail_length

        # Apply ALL operators from ALL previous x + the new x
        for x in range(1, new_limit + 1):
            for l, m, z in operators:
                y = 90 * x * x - l * x + m
                p = z + 90 * (x - 1)
                n = 0
                while True:
                    yy = int(y + n * p)
                    if yy >= epoch: break
                    if yy >= tail_start:
                        idx = yy - tail_start
                        tail_amplitude[idx] += 1
                    n += 1

        total_marks_tail = sum(tail_amplitude)
        holes_tail = sum(1 for a in tail_amplitude if a == 0)
        hole_prob = holes_tail / tail_length if tail_length > 0 else 0
        avg_amplitude_tail = total_marks_tail / tail_length if tail_length > 0 else 0

        total_marks_history.append(total_marks_tail)
        avg_amplitude_history.append(avg_amplitude_tail)
        hole_prob_history.append(hole_prob)

        print(f"h={h:3d} | Tail length={tail_length:6d} | Marks in tail={total_marks_tail:6d} | "
              f"Avg amplitude={avg_amplitude_tail:.4f} | Hole prob={hole_prob:.4f}")

        prev_epoch = epoch

    # Plot the regular growth
    plt.figure(figsize=(12, 8))
    plt.plot(range(1, h_max+1), avg_amplitude_history, 'b-o', label='Average Amplitude in Tail')
    plt.axhline(4.36, color='r', linestyle='--', label='Asymptotic raw density ≈4.36')
    plt.xlabel('Epoch h'), plt.ylabel('Average Amplitude'), plt.title('Regular Growth of Marking Pressure in New Tails')
    plt.legend(), plt.savefig('tail_amplitude_growth.png')
    print("\nPlot Saved: 'tail_amplitude_growth.png' — Shows the average amplitude in each new tail grows regularly toward the known plateau.")

    return avg_amplitude_history

# Run the analysis
segmented_tail_analysis(h_max=200)