import cmath
import math
import matplotlib.pyplot as plt

def drLD(x, l, m, z, amplitude, segment_start, segment_width):
    y = 90 * (x * x) - l * x + m
    if y < 0: y = 0
    p = z + 90 * (x - 1)
    n = 0
    while True:
        yy = int(y + n * p)
        if yy >= segment_start + segment_width: break
        if yy >= segment_start:
            idx = yy - segment_start
            amplitude[idx] += 1
        n += 1

def segmented_sieve(h_max=100, w=1000):
    densities = []
    amplitudes = []
    for h in range(1, h_max + 1):
        epoch = 90 * h * h - 12 * h + 1
        a = 90
        b = -300
        c = 250 - epoch
        d = (b**2) - (4 * a * c)
        sol2 = (-b + cmath.sqrt(d)) / (2 * a)
        new_limit = int(sol2.real) + 1

        num_segments = math.ceil(epoch / w)
        segment_hole_densities = []
        segment_total_amplitudes = []
        for s in range(num_segments):
            segment_start = s * w
            segment_end = min((s + 1) * w, int(epoch))
            segment_width = segment_end - segment_start
            amplitude = [0] * segment_width

            for x in range(1, new_limit + 1):
                drLD(x, 72, -1, 17, amplitude, segment_start, segment_width)
                drLD(x, 72, -1, 91, amplitude, segment_start, segment_width)
                drLD(x, 108, 29, 19, amplitude, segment_start, segment_width)
                drLD(x, 108, 29, 53, amplitude, segment_start, segment_width)
                drLD(x, 72, 11, 37, amplitude, segment_start, segment_width)
                drLD(x, 72, 11, 71, amplitude, segment_start, segment_width)
                drLD(x, 18, 0, 73, amplitude, segment_start, segment_width)
                drLD(x, 18, 0, 89, amplitude, segment_start, segment_width)
                drLD(x, 102, 20, 11, amplitude, segment_start, segment_width)
                drLD(x, 102, 20, 67, amplitude, segment_start, segment_width)
                drLD(x, 138, 52, 13, amplitude, segment_start, segment_width)
                drLD(x, 138, 52, 29, amplitude, segment_start, segment_width)
                drLD(x, 102, 28, 31, amplitude, segment_start, segment_width)
                drLD(x, 102, 28, 47, amplitude, segment_start, segment_width)
                drLD(x, 48, 3, 49, amplitude, segment_start, segment_width)
                drLD(x, 48, 3, 83, amplitude, segment_start, segment_width)
                drLD(x, 78, 8, 23, amplitude, segment_start, segment_width)
                drLD(x, 78, 8, 79, amplitude, segment_start, segment_width)
                drLD(x, 132, 45, 7, amplitude, segment_start, segment_width)
                drLD(x, 132, 45, 41, amplitude, segment_start, segment_width)
                drLD(x, 78, 16, 43, amplitude, segment_start, segment_width)
                drLD(x, 78, 16, 59, amplitude, segment_start, segment_width)
                drLD(x, 42, 4, 61, amplitude, segment_start, segment_width)
                drLD(x, 42, 4, 77, amplitude, segment_start, segment_width)

            total_marks = sum(amplitude)
            holes = sum(1 for amp in amplitude if amp == 0)
            density = holes / segment_width if segment_width > 0 else 0
            average_amplitude = total_marks / segment_width if segment_width > 0 else 0
            segment_hole_densities.append(density)
            segment_total_amplitudes.append(total_marks)

        # Aggregate for epoch
        avg_density = sum(segment_hole_densities) / num_segments if num_segments > 0 else 0
        total_amplitude = sum(segment_total_amplitudes)
        densities.append(avg_density)
        amplitudes.append(total_amplitude)
        print(f"h={h}, Epoch={epoch:,}, Avg Hole Density={avg_density:.4f}, Total Amplitude={total_amplitude:,}")

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, h_max + 1), amplitudes, 'b-o', label='Total Amplitude')
    plt.plot(range(1, h_max + 1), [4.36 * (90 * h * h - 12 * h + 1) for h in range(1, h_max + 1)], 'r--', label='Expected (4.36 * epoch)')
    plt.xlabel('h (Epoch)'), plt.ylabel('Total Amplitude'), plt.title('Regular Growth of Total Amplitude per Epoch')
    plt.legend(), plt.savefig('total_amplitude_growth.png')
    print("Plot Saved: 'total_amplitude_growth.png' — Shows regular quadratic growth, pressure constant.")

# Run
segmented_sieve()