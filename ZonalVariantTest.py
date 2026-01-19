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
    
    active_operators = set()  # Unique active ops
    
    intersections = []  # To record details
    
    for x in range(1, new_limit + 1):
        for op in all_ops:
            l_val, m_val, z, primitive = op
            y0 = 90 * x * x - l_val * x + m_val
            p = z + 90 * (x - 1)
            if p <= 0 or y0 >= segment_end:
                continue
            current = y0
            if y0 < segment_start:
                diff = segment_start - y0
                n = (diff + p - 1) // p
                current = y0 + n * p
            
            first_mark = None
            while current < segment_end:
                if current >= segment_start:
                    idx = current - segment_start
                    list_amp[idx] += 1
                    if first_mark is None:
                        first_mark = current
                    active_operators.add(op)
                current += p
            
            if first_mark is not None:
                intersections.append({
                    'op': op,
                    'x': x,
                    'y0': y0,
                    'p': p,
                    'first_mark': first_mark
                })
    
    total_ops = len(all_ops)
    active_ops = len(active_operators)
    avoided_ops = total_ops - active_ops
    
    # Sort intersections by first_mark
    intersections.sort(key=lambda d: d['first_mark'])
    
    return list_amp, active_ops, avoided_ops, total_ops, intersections

# Test for class [11] in span [2000000000, 2000000100]
classes = [11]
segment_start = 2000000000
segment_end = 2000000100

amp, active, avoided, total, inter = zonal_sieve_segment(classes, segment_start, segment_end)

print("Amplitudes:", amp)
print("Intersections:", inter)
print("Active:", active, "Avoided:", avoided, "Total:", total)