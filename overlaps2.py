import cmath
import math
import numpy as np
import matplotlib.pyplot as plt

def compute_marks_in_slice(h_max=100, slice_start=0, slice_width=1000):
    # Precompute new_limit for each h
    new_limits = []
    for h in range(1, h_max + 1):
        epoch = 90 * (h * h) - 12 * h + 1
        a = 90
        b = -300
        c = 250 - epoch
        d = (b**2) - (4 * a * c)
        sol2 = (-b + cmath.sqrt(d)) / (2 * a)
        new_limits.append(int(sol2.real) + 1)

    # Cumulative marks in slice (amplitude sum 0-999)
    marks_history = []
    amplitude = np.zeros(slice_width, dtype=int)  # Only track 0-999

    def drLD(x, l, m, z):
        y = 90 * (x * x) - l * x + m
        if y >= slice_width: return  # Beyond slice
        p = z + 90 * (x - 1)
        n = 0
        while True:
            yy = int(y + n * p)
            if yy >= slice_width: break
            if yy >= 0:
                amplitude[yy] += 1
            n += 1

    for h in range(1, h_max + 1):
        x_start = new_limits[h-2] if h > 1 else 0  # New x for this h
        x_end = new_limits[h-1]
        for x in range(x_start + 1, x_end + 1):
            # Apply all 24 operators for new x
            drLD(x, 72, -1, 17)
            drLD(x, 72, -1, 91)
            drLD(x, 108, 29, 19)
            drLD(x, 108, 29, 53)
            drLD(x, 72, 11, 37)
            drLD(x, 72, 11, 71)
            drLD(x, 18, 0, 73)
            drLD(x, 18, 0, 89)
            drLD(x, 102, 20, 11)
            drLD(x, 102, 20, 67)
            drLD(x, 138, 52, 13)
            drLD(x, 138, 52, 29)
            drLD(x, 102, 28, 31)
            drLD(x, 102, 28, 47)
            drLD(x, 48, 3, 49)
            drLD(x, 48, 3, 83)
            drLD(x, 78, 8, 23)
            drLD(x, 78, 8, 79)
            drLD(x, 132, 45, 7)
            drLD(x, 132, 45, 41)
            drLD(x, 78, 16, 43)
            drLD(x, 78, 16, 59)
            drLD(x, 42, 4, 61)
            drLD(x, 42, 4, 77)
        total_marks = np.sum(amplitude)
        marks_history.append(total_marks)
        print(f"h={h}, Total Marks in 0-999: {total_marks}, Density: {total_marks / slice_width:.2f}")

    # Plot Growth
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, h_max + 1), marks_history, 'b-o', label='Total Marks')
    plt.axhline(4360, color='r', linestyle='--', label='Expected Plateau (~4.36 * 1000)')
    plt.xlabel('h (Epoch)'), plt.ylabel('Total Marks in 0-999'), plt.title('Marks Growth in Fixed Slice with Epoch')
    plt.legend(), plt.savefig('marks_growth.png')
    print("Plot Saved: 'marks_growth.png' — Shows growth to plateau, pressure bounded.")
    return marks_history

# Run
compute_marks_in_slice()