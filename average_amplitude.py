import cmath
import math
import matplotlib.pyplot as plt

def compute_average_amplitude(h_max=1000):
    averages = []
    for h in range(1, h_max + 1):
        epoch = 90 * h * h - 12 * h + 1
        a = 90
        b = -300
        c = 250 - epoch
        d = (b**2) - (4 * a * c)
        sol2 = (-b + cmath.sqrt(d)) / (2 * a)
        new_limit = int(sol2.real) + 1

        amplitude = [0] * int(epoch + 100)

        def drLD(x, l, m, z, listvar):
            y = 90 * (x * x) - l * x + m
            if y < 0 or y >= len(listvar): return
            listvar[int(y)] += 1
            p = z + 90 * (x - 1)
            n = 1
            while True:
                yy = int(y + n * p)
                if yy >= len(listvar): break
                listvar[yy] += 1
                n += 1

        for x in range(1, new_limit + 1):
            drLD(x, 72, -1, 17, amplitude)
            drLD(x, 72, -1, 91, amplitude)
            drLD(x, 108, 29, 19, amplitude)
            drLD(x, 108, 29, 53, amplitude)
            drLD(x, 72, 11, 37, amplitude)
            drLD(x, 72, 11, 71, amplitude)
            drLD(x, 18, 0, 73, amplitude)
            drLD(x, 18, 0, 89, amplitude)
            drLD(x, 102, 20, 11, amplitude)
            drLD(x, 102, 20, 67, amplitude)
            drLD(x, 138, 52, 13, amplitude)
            drLD(x, 138, 52, 29, amplitude)
            drLD(x, 102, 28, 31, amplitude)
            drLD(x, 102, 28, 47, amplitude)
            drLD(x, 48, 3, 49, amplitude)
            drLD(x, 48, 3, 83, amplitude)
            drLD(x, 78, 8, 23, amplitude)
            drLD(x, 78, 8, 79, amplitude)
            drLD(x, 132, 45, 7, amplitude)
            drLD(x, 132, 45, 41, amplitude)
            drLD(x, 78, 16, 43, amplitude)
            drLD(x, 78, 16, 59, amplitude)
            drLD(x, 42, 4, 61, amplitude)
            drLD(x, 42, 4, 77, amplitude)

        amplitude = amplitude[:int(epoch)]
        total_marks = sum(amplitude)
        average = total_marks / epoch
        averages.append(average)
        print(f"h={h}, Epoch={epoch:,}, Total Marks={total_marks:,}, Average Amplitude={average:.4f}")

    # Plot Growth
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, h_max + 1), averages, 'b-o', label='Average Amplitude')
    plt.axhline(4.36, color='r', linestyle='--', label='Asymptotic Raw Density (~4.36)')
    plt.xlabel('h (Epoch)'), plt.ylabel('Average Amplitude'), plt.title('Regular Growth of Average Amplitude per Epoch')
    plt.legend(), plt.savefig('average_amplitude_growth.png')
    print("Plot Saved: 'average_amplitude_growth.png' — Shows regular growth to plateau, pressure stabilizing.")

# Run
compute_average_amplitude()