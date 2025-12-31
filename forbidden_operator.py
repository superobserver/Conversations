import math
import cmath

def check_overlap(l, m, z, new_limit, start_y, end_y):
    for x in range(1, int(new_limit) + 1):
        y = 90 * (x * x) - l * x + m
        p = z + (90 * (x - 1))
        if p <= 0:
            continue
        # Check if AP intersects span
        if y >= start_y and y <= end_y:
            return True
        min_n = math.ceil((start_y - y) / p) if y < start_y else 0
        first = y + min_n * p
        if first <= end_y:
            return True
    return False

def pre_tabulate_forbidden(h=50, start_y=1000, end_y=1999):
    epoch = 90 * (h * h) - 12 * h + 1
    a = 90
    b = -300
    c = 250 - epoch
    d = (b ** 2) - (4 * a * c)
    sol2 = (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = sol2.real

    operators = [
        (72, -1, 17), (72, -1, 91),
        (108, 29, 19), (108, 29, 53),
        (72, 11, 37), (72, 11, 71),
        (18, 0, 73), (18, 0, 89),
        (102, 20, 11), (102, 20, 67),
        (138, 52, 13), (138, 52, 29),
        (102, 28, 31), (102, 28, 47),
        (48, 3, 49), (48, 3, 83),
        (78, 8, 23), (78, 8, 79),
        (132, 45, 7), (132, 45, 41),
        (78, 16, 43), (78, 16, 59),
        (42, 4, 61), (42, 4, 77)
    ]

    forbidden = []
    for l, m, z in operators:
        if not check_overlap(l, m, z, new_limit, start_y, end_y):
            forbidden.append((l, m, z))

    print("Forbidden operators (no overlap):", len(forbidden))
    print("Sample forbidden:", forbidden[:3])

    # Optional: Full marking to verify no marks from forbidden
    full_amp = [0] * (end_y - start_y + 1)
    for l, m, z in operators:
        for x in range(1, int(new_limit) + 1):
            y = 90 * (x * x) - l * x + m
            p = z + (90 * (x - 1))
            if p <= 0:
                continue
            n = 0
            while True:
                next_y = y + n * p
                if next_y > end_y:
                    break
                if start_y <= next_y <= end_y:
                    full_amp[next_y - start_y] += 1
                n += 1

    # Optimized: Mark only non-forbidden
    opt_amp = [0] * (end_y - start_y + 1)
    active = [op for op in operators if op not in forbidden]
    for l, m, z in active:
        for x in range(1, int(new_limit) + 1):
            y = 90 * (x * x) - l * x + m
            p = z + (90 * (x - 1))
            if p <= 0:
                continue
            n = 0
            while True:
                next_y = y + n * p
                if next_y > end_y:
                    break
                if start_y <= next_y <= end_y:
                    opt_amp[next_y - start_y] += 1
                n += 1

    match = full_amp == opt_amp
    print("Amplitudes match:", match)
    print("Active operators:", len(active))

if __name__ == "__main__":
    pre_tabulate_forbidden()
