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
    if segment_start >= segment_end:
        raise ValueError("segment_start must be less than segment_end")
    
    a = 90
    b = -300
    c = 250 - segment_end
    d = (b**2) - (4 * a * c)
    sol2 = (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = int(sol2.real) + 1
    
    segment_width = segment_end - segment_start
    list_amp = [0] * segment_width
    
    all_ops = []
    for k in classes:
        all_ops.extend(get_operators(k))
    
    active_ops = 0
    
    for x in range(1, new_limit + 1):
        for l_val, m_val, z, primitive in all_ops:
            y0 = 90 * x * x - l_val * x + m_val
            p = z + 90 * (x - 1)
            if p <= 0:
                continue
            
            # CRT-based pruning: Check if AP y0 mod p intersects [start, end)
            # First possible n >= start: current = y0 + ceil((start - y0)/p) * p
            diff = segment_start - y0
            n_start = (diff + p - 1) // p if diff > 0 else 0
            first_n = y0 + n_start * p
            if first_n >= segment_end:
                continue  # No intersection: prune
            
            # Last possible n < end
            last_n = y0 + ((segment_end - 1 - y0) // p) * p
            if first_n > last_n:
                continue  # Empty range: prune
            
            # Last digit binding (mod 10 prune): If class k last digit fixed, but for composite mark,
            # optional: check if y0 mod 10 allows composite last digit, but skip for now as loose
            
            # If passes, proceed to mark (no while if no intersection, but here we know there is)
            current = first_n
            marked = False
            while current < segment_end:
                idx = current - segment_start
                list_amp[idx] += 1
                marked = True
                current += p
            
            if marked:
                active_ops += 1  # But this overcounts; fix to unique in full theory
    
    total_ops = len(all_ops)
    avoided_ops = total_ops - active_ops
    return list_amp, active_ops, avoided_ops, total_ops


"""

def zonal_sieve_segment(classes, segment_start, segment_end):
    ""
    Modified zonal sieve: Mark only a single segment [start, end), return amplitude list for it.
    Processes operators intersecting the segment, for local analysis (e.g., adjacencies).
    Handles large indices with integer arithmetic.
    Now reports active operators (those that mark at least once in the span).
    ""
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
    
    active_ops = 0  # Count of operators that actually mark in the span
    
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
            
            marked = False  # Flag for this operator at this x
            while current < segment_end:
                if current >= segment_start:
                    idx = current - segment_start
                    list_amp[idx] += 1
                    marked = True
                current += p
            
            if marked:
                active_ops += 1  # Increment only if this op marked at least once
    
    total_ops = len(all_ops)
    avoided_ops = total_ops - active_ops
    return list_amp, active_ops, avoided_ops, total_ops
"""
def process_class(args):
    k, segment_start, segment_end, oeis_id, description = args
    class_tag = f"[classes = [{k}]]"
    classes = [k]
    print(f"{class_tag} {oeis_id} {description}")
    print(f"{class_tag} Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
    amp_segment, active_ops, avoided_ops, total_ops = zonal_sieve_segment(classes, segment_start, segment_end)
    print(f"{class_tag} Segment amplitudes (first 100): {amp_segment[:100]}")
    oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
    print(f"{class_tag} OEIS-like sequence (hole indices): {oeis_sequence}")
    holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
    print(f"{class_tag} Holes in segment: {holes_in_segment}")
    print(f"{class_tag} Operator report: Active = {active_ops}, Avoided = {avoided_ops}, Total = {total_ops}\n")
    # Return for potential collection
    return (k, holes_in_segment, oeis_sequence, active_ops, avoided_ops)

if __name__ == "__main__":
    segment_start = 2000000000000000
    segment_end =   2000000000001000
    
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
    
    # Aggregate total holes and operator stats across all classes (for global insights)
    total_holes = sum(res[1] for res in results)
    total_active = sum(res[3] for res in results)
    total_avoided = sum(res[4] for res in results)
    print(f"Global report: Total holes across all classes: {total_holes}")
    print(f"Total active operators: {total_active}, Total avoided: {total_avoided}")