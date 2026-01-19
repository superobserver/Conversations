import cmath
import multiprocessing

def get_operators(k):
    R = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59, 61, 67, 71, 73, 77, 79, 83, 89]
    operators = []
    seen_pairs = set()
    for z in R:
        try:
            o = (k * pow(z, -1, 90)) % 90
            if o not in R:
                continue
            pair = tuple(sorted([z, o]))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            z_eff = 91 if z == 1 else z
            o_eff = 91 if o == 1 else o
            l = 180 - (z_eff + o_eff)
            m = 90 - (z_eff + o_eff) + (z_eff * o_eff - k) // 90
            operators.append((l, m, z_eff, k))
            if z != o:
                operators.append((l, m, o_eff, k))
        except ValueError:
            continue
    return operators

def zonal_sieve_segment(classes, segment_start, segment_end):
    """
    Modified zonal sieve: Mark only a single segment [start, end), return amplitude list for it.
    Processes operators intersecting the segment, for local analysis (e.g., adjacencies).
    Handles large indices with integer arithmetic.
    """
    if segment_start >= segment_end:
        raise ValueError("segment_start must be less than segment_end")
    
    a = 90
    b = -300
    c = 250 - segment_end  # Bound based on segment_end for pruning x iterations
    d = (b**2) - (4 * a * c)
    sol2 = (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = int(sol2.real) + 1
    
    segment_width = segment_end - segment_start
    list_amp = [0] * segment_width  # Amplitudes for the segment
    
    all_ops = []
    for k in classes:
        all_ops.extend(get_operators(k))
    
    for x in range(1, new_limit + 1):
        for l_val, m_val, z, primitive in all_ops:
            y0 = 90 * x * x - l_val * x + m_val
            p = z + 90 * (x - 1)
            if p <= 0 or y0 >= segment_end:
                continue
            if y0 < segment_start:
                diff = segment_start - y0
                # Integer ceiling division to handle large numbers precisely
                n = (diff + p - 1) // p
                current = y0 + n * p
            else:
                current = y0
            while current < segment_end:
                if current >= segment_start:
                    idx = current - segment_start
                    list_amp[idx] += 1
                current += p
    return list_amp

def process_class(args):
    k, segment_start, segment_end, oeis_id, description = args
    classes = [k]
    print(f"{oeis_id} {description}")
    print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
    amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
    print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
    oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
    print("OEIS-like sequence (hole indices):", oeis_sequence)
    holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
    print(f"Holes in segment: {holes_in_segment}\n")
    # Return for potential collection, though we print directly
    return (k, holes_in_segment, oeis_sequence)

if __name__ == "__main__":
    segment_start = 10000
    segment_end =   11000
    
    # List of all classes with their OEIS IDs and descriptions
    class_list = [
        (7, "A202110", "Numbers k such that 90*k + 7 is prime."),
        (11, "A201804", "Numbers k such that 90*k + 11 is prime."),
        (13, "A201816", "Numbers k such that 90*k + 13 is prime."),
        (17, "A202115", "Numbers k such that 90*k + 17 is prime."),
        (19, "A202121", "Numbers k such that 90*k + 19 is prime."),
        (23, "A201823", "Numbers k such that 90*k + 23 is prime."),
        (29, "A201806", "Numbers k such that 90*k + 29 is prime."),
        (31, "A201810", "Numbers k such that 90*k + 31 is prime."),
        (37, "A201829", "Numbers k such that 90*k + 37 is prime."),
        (41, "A201807", "Numbers k such that 90*k + 41 is prime."),
        (43, "A201811", "Numbers k such that 90*k + 43 is prime."),
        (47, "A202124", "Numbers k such that 90*k + 47 is prime."),
        (49, "A202120", "Numbers k such that 90*k + 49 is prime."),
        (53, "A201824", "Numbers k such that 90*k + 53 is prime."),
        (59, "A201808", "Numbers k such that 90*k + 59 is prime."),
        (61, "A201812", "Numbers k such that 90*k + 61 is prime."),
        (67, "A201817", "Numbers k such that 90*k + 67 is prime."),
        (71, "A202129", "Numbers n such that 90n + 71 is prime."),
        (73, "A195993", "Numbers n such that 90n + 73 is prime."),
        (77, "A201822", "Numbers k such that 90*k + 77 is prime."),
        (79, "A202112", "Numbers k such that 90*k + 79 is prime."),
        (83, "A196007", "Numbers n such that 90n + 83 is prime."),
        (89, "A202116", "Numbers k such that 90*k + 89 is prime."),
        (91, "A224889", "Numbers n such that 90n + 91 is prime."),  # Note: 91 mod 90 = 1
        (1, "A181732", "Numbers n such that 90n + 1 is prime.")
    ]
    
    # Prepare arguments for each class
    args_list = [(k, segment_start, segment_end, oeis_id, desc) for k, oeis_id, desc in class_list]
    
    # Use multiprocessing Pool for parallel execution
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = pool.map(process_class, args_list)
    
    # Optionally process results (e.g., aggregate hole counts), but since we print in processes, this is for completeness
    total_holes = sum(res[1] for res in results)
    print(f"Total holes across all classes: {total_holes}")