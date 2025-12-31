import math
import cmath

def drLD(x, l, m, z, listvar):
    y = 90 * (x * x) - l * x + m
    if y >= len(listvar) or y < 0:
        return
    listvar[int(y)] += 1
    p = z + (90 * (x - 1))
    if p == 0:
        return  # Avoid division by zero, though unlikely
    for n in range(1, int(((len(listvar) - 1 - y) / p) + 1)):
        next_y = y + (p * n)
        if next_y >= len(listvar):
            break
        listvar[int(next_y)] += 1

def sieve_with_hilbert(h=10, residue_k=17):
    epoch = 90 * (h * h) - 12 * h + 1
    a = 90
    b = -300
    c = 250 - epoch
    d = (b ** 2) - (4 * a * c)
    sol2 = (-b + cmath.sqrt(d)) / (2 * a)
    new_limit = sol2.real  # Use real part

    list_amp = [0] * (epoch + 100)

    # 24 operators for class 17 (from the manuscript/code)
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
            drLD(x, l, m, z, list_amp)

    list_amp = list_amp[:epoch]  # Trim to exact epoch

    # Hilbert chains: list of (y, 1/p, skew_phase) per operator
    chains = [[] for _ in operators]
    for op_idx, (l, m, z) in enumerate(operators):
        for x in range(1, int(new_limit) + 1):
            y = 90 * (x * x) - l * x + m
            if 0 <= y < epoch:
                p = z + 90 * (x - 1)
                if p != 0:
                    one_p = 1.0 / p
                    skew_phase = (l * x + m) / 90.0  # Normalized skew
                    chains[op_idx].append((y, one_p, skew_phase))

    # Wavefunction approximation (real part sum)
    wave = [0.0] * epoch
    for chain in chains:
        for y, one_p, skew_phase in chain:
            wave[int(y)] += one_p * math.cos(2 * math.pi * skew_phase)

    # State tuples with observables and confidence
    states = []
    for y in range(epoch):
        n = 90 * y + residue_k
        amp = list_amp[y]
        is_prime = (amp == 0)
        last_digit = n % 10
        digital_root = n % 9 if n % 9 != 0 else 9
        confidence = 1.0 if is_prime else 1.0 - 1.0 / (amp + 1)  # Range space knowledge
        states.append({
            'y_index': y,
            'number_n': n,
            'amplitude': amp,
            'prime_property': is_prime,
            'last_digit_observable': last_digit,
            'digital_root_observable': digital_root,
            'range_space_confidence': confidence,
            'wave_value': wave[y]
        })

    return states, chains, wave

# Example run and test
if __name__ == "__main__":
    states, chains, wave = sieve_with_hilbert(h=10)
    print("Total holes (primes):", sum(1 for s in states if s['prime_property']))
    print("Sample states (first 10):")
    for s in states[:10]:
        print(s)
    print("\nSample chain 0 (first 5 entries):", chains[0][:5])
    print("Sample wave values (first 10 y):", wave[:10])