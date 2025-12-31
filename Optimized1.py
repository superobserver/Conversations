import math
import cmath

def drLD(x, l, m, z, listvar, start_y=0, end_y=None):
    y = 90 * (x * x) - l * x + m
    if end_y is None:
        end_y = len(listvar) - 1
    p = z + (90 * (x - 1))
    if p <= 0:
        return 0
    # Skip if no intersection with span
    if y > end_y:
        return 0  # y too large, no forward marks
    min_n = 0
    if y < start_y:
        min_n = math.ceil((start_y - y) / p)
    first_mark = y + min_n * p
    if first_mark > end_y:
        return 0  # Cannot reach span
    marks_done = 0
    n = min_n
    while True:
        next_y = y + n * p
        if next_y > end_y:
            break
        if 0 <= next_y < len(listvar):
            listvar[int(next_y)] += 1
            marks_done += 1
        n += 1
    return marks_done

def quadratic_sieve_span_optimized(h=50, residue_k=17, start_y=1000, end_y=1999):
    epoch = 90 * (h * h) - 12 * h + 1
    a = 90
    b = -300
    c = 250 - epoch
    d = (b ** 2) - (4 * a * c)
    sol2 = (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = sol2.real

    # Full dataset for comparison (original NewGrok17.py style)
    full_amp = [0] * int(epoch + 100)
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

    for x in range(1, int(new_limit) + 1):
        for l, m, z in operators:
            drLD(x, l, m, z, full_amp, 0, epoch - 1)

    full_amp = full_amp[:epoch]

    # Optimized for span: Use same drLD but with span bounds, count marks
    span_amp = [0] * (end_y - start_y + 1)
    total_marks = 0
    skipped_calls = 0
    for x in range(1, int(new_limit) + 1):
        for l, m, z in operators:
            marks = drLD(x, l, m, z, span_amp, start_y, end_y)
            if marks == 0:
                skipped_calls += 1
            total_marks += marks

    # Test match: Extract full span for comparison
    full_span_amp = full_amp[start_y:end_y + 1]
    match = full_span_amp == span_amp

    full_holes = sum(1 for amp in full_span_amp if amp == 0)
    opt_holes = sum(1 for amp in span_amp if amp == 0)

    print("Full holes in span:", full_holes)
    print("Optimized holes in span:", opt_holes)
    print("Amplitudes match:", match)
    print("Total marks in optimized:", total_marks)
    print("Skipped drLD calls:", skipped_calls)

if __name__ == "__main__":
    quadratic_sieve_span_optimized()
