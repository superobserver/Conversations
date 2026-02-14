import math
import numpy as np
from midiutil import MIDIFile
import argparse
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# Constants (unchanged)
ALL_CLASSES = [1,7,11,13,17,19,23,29,31,37,41,43,47,49,53,59,61,67,71,73,77,79,83,89]

TWIN_PAIRS = [
    (11, 13), (17, 19), (29, 31), (41, 43), (47, 49),
    (59, 61), (71, 73), (77, 79), (89, 91)
]

def get_operators(k):
    operators = []
    seen = set()
    for z in ALL_CLASSES:
        try:
            o = (k * pow(z, -1, 90)) % 90
            if o not in ALL_CLASSES: continue
            pair = tuple(sorted([z, o]))
            if pair in seen: continue
            seen.add(pair)
            z_eff = 91 if z == 1 else z
            o_eff = 91 if o == 1 else o
            l = 180 - (z_eff + o_eff)
            m = 90 - (z_eff + o_eff) + (z_eff * o_eff - k) // 90
            operators.append((l, m, z_eff, z))
            if z != o:
                operators.append((l, m, o_eff, o))
        except ValueError:
            continue
    return operators


####### new function #############

def active_xs_for_operator(l, m, z_eff, start, end):
    """
    Conservative but correct: return x_min=1 to a safe upper bound,
    then rely on per-x hit check to filter.
    This ensures we never miss small-x generators when start is large.
    """
    hi = int(math.sqrt(end / 90)) + 500  # Extra margin for safety
    # Also allow larger x if p_x is small enough to hit multiple times
    hi = max(hi, (start // 90) + 1000)  # Rough scale for large start

    # Minimal x is always 1 (we will check if it hits)
    x_min = 1
    x_max = hi

    # Quick reject if even largest x has y(x) way above end
    y_max = 90 * x_max * x_max - l * x_max + m
    if y_max < start:
        return None

    return x_min, x_max


########## close new function ############


"""
def active_xs_for_operator(l, m, z_eff, start, end):
    lo, hi = 1, int(math.sqrt(max(end, 100000) / 90)) + 100  # Broader for high start
    # Largest x where y(x) < end
    x_max_bound = None
    left, right = lo, hi
    while left <= right:
        mid = (left + right) // 2
        y = 90 * mid * mid - l * mid + m
        if y < end:
            x_max_bound = mid
            left = mid + 1
        else:
            right = mid - 1
    if x_max_bound is None:
        return None

    # Smallest x where hit < end
    left, right = 1, x_max_bound
    x_min_hit = None
    while left <= right:
        mid = (left + right) // 2
        y = 90 * mid * mid - l * mid + m
        p = z_eff + 90 * (mid - 1)
        if p <= 0:
            left = mid + 1
            continue
        ceil_val = math.ceil((start - y) / p) if (start - y) > 0 else 0
        hit = y + p * ceil_val
        if hit < end:
            x_min_hit = mid
            right = mid - 1
        else:
            left = mid + 1

    if x_min_hit is None:
        return None

    x_min = max(1, x_min_hit)
    x_max = x_max_bound

    if x_min > x_max:
        return None

    return x_min, x_max
"""
def pre_compute_cohort(k, start, end):
    cohort = []
    ops = get_operators(k)
    for l, m, z_eff, z in ops:
        bounds = active_xs_for_operator(l, m, z_eff, start, end)
        if bounds is None:
            continue
        x_min, x_max = bounds
        for x in range(x_min, x_max + 1):
            y0 = 90 * x * x - l * x + m
            p = z_eff + 90 * (x - 1)
            if p <= 0: continue
            ceil_val = math.ceil((start - y0) / p) if (start - y0) > 0 else 0
            hit = y0 + p * ceil_val
            if start <= hit < end:
                cohort.append((k, z, x, y0, p))
    return cohort

def compute_amplitude(k, start, end):
    cohort = pre_compute_cohort(k, start, end)
    print(f"  Class {k:2d}: cohort size = {len(cohort)}")
    ###### new addition ######
    print(f"  Class {k:2d}: cohort size = {len(cohort)}")
    if len(cohort) == 0:
        print(f"  WARNING: No markings possible in [{start}, {end}) for class {k}")
    ###### end new additon ######
    max_amp = max(amp) if amp else 0
    print(f"  Max amplitude = {max_amp}, zeros = {sum(1 for a in amp if a == 0)}")
    if not cohort:
        print(f"  WARNING: No active operators for class {k} in range -- all unmarked")
    width = end - start
    amp = np.zeros(width, dtype=np.int32)  # Higher dtype for dense marks
    for _, _, _, y0, p in cohort:
        current = y0
        while current < start:
            current += p
        while current < end:
            idx = current - start
            if 0 <= idx < width:
                amp[idx] += 1
            current += p
    return amp.tolist()

# Parallel worker (unchanged)
def worker_compute_amp(args):
    k, start, end = args
    amp = compute_amplitude(k, start, end)
    return k, amp

# Main (unchanged except segmented flag integration)
def generate_midi(mode='twin', start=0, end=100000, time_scale=0.008,
                  note_duration=0.25, base_pitch=60, pitch_sep=7,
                  output_file="sieve_music.mid", parallel=True, output_amp_lists=False, segmented=False):

    if segmented:
        # In segmented mode, broaden hi in active_xs if needed
        pass  # Already broadened in function

    # Rest of the function unchanged...

# CLI (add --segmented flag)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Quadratic Sieve MIDI Generator (parallel & pruned)")
    parser.add_argument('--mode', choices=['twin', 'full', 'composite'], default='twin')
    parser.add_argument('--start', type=int, default=10000)
    parser.add_argument('--end', type=int, default=100000)
    parser.add_argument('--time_scale', type=float, default=0.008)
    parser.add_argument('--note_duration', type=float, default=0.25)
    parser.add_argument('--base_pitch', type=int, default=60)
    parser.add_argument('--pitch_sep', type=int, default=7)
    parser.add_argument('--output', type=str, default="sieve_music99999.mid")
    parser.add_argument('--no-parallel', action='store_true')
    parser.add_argument('--output_amp', action='store_true')
    parser.add_argument('--segmented', action='store_true', help="Treat as high segment, force broader x search")

    args = parser.parse_args()
    generate_midi(
        mode=args.mode,
        start=args.start,
        end=args.end,
        time_scale=args.time_scale,
        note_duration=args.note_duration,
        base_pitch=args.base_pitch,
        pitch_sep=args.pitch_sep,
        output_file=args.output,
        parallel=not args.no_parallel,
        output_amp_lists=args.output_amp,
        segmented=args.segmented
    )