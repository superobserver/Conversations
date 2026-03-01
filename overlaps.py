import cmath
import math

def compute_overlapping(h_max=100, range_size=1000):
    results = []
    for h in range(1, h_max + 1):
        epoch = 90 * (h * h) - 12 * h + 1
        a = 90
        b = -300
        c = 250 - epoch
        d = (b**2) - (4 * a * c)
        sol2 = (-b + cmath.sqrt(d)) / (2 * a)
        new_limit = int(sol2.real) + 1

        overlapping = set()
        def drLD(x, l, m, z):
            y = 90 * (x * x) - l * x + m
            if y >= range_size:
                return
            p = z + 90 * (x - 1)
            current = y
            while current < range_size:
                if 0 <= current < range_size:
                    overlapping.add((z, x))
                current += p

        for x in range(1, new_limit):
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

        num_overlapping = len(overlapping)
        results.append((h, num_overlapping))
    return results

overlaps = compute_overlapping()
print(overlaps)