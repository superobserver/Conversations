import cmath

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

# Example: Sieve a small range at a large index for class [17]
segment_start = 1000000
segment_end =   1001000




classes = [7] 
print("A202110 Numbers k such that 90*k + 7 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")


classes = [11] 
print("A201804 Numbers k such that 90*k + 11 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")


classes = [13] 
print("A201816 Numbers k such that 90*k + 13 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [17] 
print("A202115 Numbers k such that 90*k + 17 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [19] 
print("A196000 Numbers k such that 90*k + 19 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [23] 
print("A201820 Numbers k such that 90*k + 23 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [29] 
print("A201739 Numbers n such that 90*n + 29 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [31] 
print("A201819 Numbers n such that 90*n + 31 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [37] 
print("A198382 Numbers n such that 90n + 37 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [41] 
print("A202104 Numbers k such that 90*k + 41 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [43] 
print("A202105 Numbers k such that 90*k + 43 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")
classes = [47] 
print("A201734 Numbers n such that 90*n + 47 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [49] 
print("A201818 Numbers k such that 90*k + 49 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [53] 
print("A202114 Numbers k such that 90*k + 53 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [59] 
print("A202101 Numbers k such that 90*k + 59 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [61] 
print("A202113 Numbers k such that 90*k + 61 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [67] 
print("A201817 Numbers k such that 90*k + 67 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [71] 
print("A202129 Numbers n such that 90n + 71 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [73] 
print("A195993 Numbers n such that 90n + 73 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [77] 
print("A201822 Numbers k such that 90*k + 77 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [79] 
print("A202112 Numbers k such that 90*k + 79 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [83] 
print("A196007 Numbers n such that 90n + 83 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [89] 
print("A202116 Numbers k such that 90*k + 89 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [91] 
print("A224889 Numbers n such that 90n + 91 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")

classes = [1] 
print("A181732 Numbers n such that 90n + 1 is prime.")
print(f"Sieving segment [{segment_start}, {segment_end}) for classes {classes}")
amp_segment = zonal_sieve_segment(classes, segment_start, segment_end)
print(f"Segment amplitudes (first 100): {amp_segment[:100]}")
oeis_sequence = [segment_start + i for i, x in enumerate(amp_segment) if x == 0]
print("OEIS-like sequence (hole indices):", oeis_sequence)
holes_in_segment = sum(1 for amp in amp_segment if amp == 0)
print(f"Holes in segment: {holes_in_segment}")