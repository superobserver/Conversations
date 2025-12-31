import math
import cmath

def drLD(x, l, m, z, listvar):
    y = 90 * (x * x) - l * x + m
    if y >= len(listvar) or y < 0:
        return
    listvar[int(y)] += 1
    p = z + (90 * (x - 1))
    if p == 0:
        return
    for n in range(1, int(((len(listvar) - 1 - y) / p) + 1)):
        next_y = y + (p * n)
        if next_y >= len(listvar):
            break
        listvar[int(next_y)] += 1

def drLD_optimized(x, l, m, z, span_amp, start_y, end_y):
    y = 90 * (x * x) - l * x + m
    p = z + (90 * (x - 1))
    if p == 0:
        return 0
    min_n = 0 if y >= start_y else math.ceil((start_y - y) / p)
    first_mark = y + min_n * p
    if first_mark > end_y:
        return 0
    marks = 0
    n = min_n
    while True:
        next_y = y + n * p
        if next_y > end_y:
            break
        local_idx = next_y - start_y
        if 0 <= local_idx < len(span_amp):
            span_amp[local_idx] += 1
            marks += 1
        n += 1
    return marks

def sieve_full_and_optimized_span(h=50, residue_k=17, start_y=1000, end_y=1099):
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

    # Full sieve (original NewGrok17.py)
    full_amp = [0] * int(epoch + 100)
    full_marks = 0
    for x in range(1, int(new_limit) + 1):
        for l, m, z in operators:
            drLD(x, l, m, z, full_amp)
            full_marks += 1  # Count calls (operations proxy)

    full_amp = full_amp[:epoch]
    full_span_amp = full_amp[start_y:end_y + 1]
    full_holes = sum(1 for amp in full_span_amp if amp == 0)

    # Optimized span sieve
    span_amp = [0] * (end_y - start_y + 1)
    opt_marks = 0
    for x in range(1, int(new_limit) + 1):
        for l, m, z in operators:
            marks = drLD_optimized(x, l, m, z, span_amp, start_y, end_y)
            opt_marks += marks if marks > 0 else 0  # Count effective marks

    opt_holes = sum(1 for amp in span_amp if amp == 0)

    match = full_span_amp == span_amp

    print("Full holes in span:", full_holes)
    print("Optimized holes in span:", opt_holes)
    print("Amplitudes match:", match)
    print("Full operations (calls):", full_marks)
    print("Optimized operations (marks):", opt_marks)

if __name__ == "__main__":
    sieve_full_and_optimized_span()
