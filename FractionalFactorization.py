import math
import numpy as np
from fractions import Fraction
from sympy import continued_fraction

# Operator from papers
l = 102
m = 20
z_eff = 11
k = 11

def y_x(x):
    return 90 * x**2 - l * x + m

def p_x(x):
    return z_eff + 90 * (x - 1)

def f_x(n, x):
    y = y_x(x)
    p = p_x(x)
    if p <= 0:
        return None
    raw_b = (n - y) / p
    return raw_b - math.floor(raw_b)

def max_x_for_n(n):
    return int(math.sqrt((90 * n + k - m + 1) / 90)) + 1

# Hole n=3 (281 prime)
print("Hole n=3:")
max_x_h = max_x_for_n(3)
f_h = []
rates_h = []
cf_lengths_h = []
for x in range(1, max_x_h + 1):
    f = f_x(3, x)
    if f is not None:
        f_h.append(f)
        if len(f_h) > 1:
            rates_h.append(f_h[-1] - f_h[-2])
        rat = Fraction.from_float(f).limit_denominator(1000)
        cf = continued_fraction(rat)
        cf_lengths_h.append(len(list(cf)))
print("f:", f_h)
print("rates:", rates_h)
print("CF len:", cf_lengths_h)
print("rate var:", np.var(rates_h) if rates_h else "N/A")

# Chained n=4 (371 composite)
print("\nChained n=4:")
max_x_c = max_x_for_n(4)
f_c = []
rates_c = []
cf_lengths_c = []
for x in range(1, max_x_c + 1):
    f = f_x(4, x)
    if f is not None:
        f_c.append(f)
        if len(f_c) > 1:
            rates_c.append(f_c[-1] - f_c[-2])
        rat = Fraction.from_float(f).limit_denominator(1000)
        cf = continued_fraction(rat)
        cf_lengths_c.append(len(list(cf)))
print("f:", f_c)
print("rates:", rates_c)
print("CF len:", cf_lengths_c)
print("rate var:", np.var(rates_c) if rates_c else "N/A")

# Larger hole n=20 (1811 prime)
print("\nLarge hole n=20:")
max_x_lh = max_x_for_n(20)
f_lh = []
rates_lh = []
cf_lengths_lh = []
for x in range(1, max_x_lh + 1):
    f = f_x(20, x)
    if f is not None:
        f_lh.append(f)
        if len(f_lh) > 1:
            rates_lh.append(f_lh[-1] - f_lh[-2])
        rat = Fraction.from_float(f).limit_denominator(1000)
        cf = continued_fraction(rat)
        cf_lengths_lh.append(len(list(cf)))
print("f:", f_lh)
print("rates:", rates_lh)
print("CF len:", cf_lengths_lh)
print("rate var:", np.var(rates_lh) if rates_lh else "N/A")

# Larger chained n=22 (1991 composite)
print("\nLarge chained n=22:")
max_x_lc = max_x_for_n(22)
f_lc = []
rates_lc = []
cf_lengths_lc = []
for x in range(1, max_x_lc + 1):
    f = f_x(22, x)
    if f is not None:
        f_lc.append(f)
        if len(f_lc) > 1:
            rates_lc.append(f_lc[-1] - f_lc[-2])
        rat = Fraction.from_float(f).limit_denominator(1000)
        cf = continued_fraction(rat)
        cf_lengths_lc.append(len(list(cf)))
print("f:", f_lc)
print("rates:", rates_lc)
print("CF len:", cf_lengths_lc)
print("rate var:", np.var(rates_lc) if rates_lc else "N/A")
