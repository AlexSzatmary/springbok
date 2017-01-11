from collections import defaultdict
import glob
import numpy as np


def find():
    d_L = defaultdict(lambda: [])
    for f in glob.glob('*-range.txt'):
        a = np.loadtxt(f)
        (gamma_L, phi_E, k, r_L, range_) = a
        d_L[gamma_L, phi_E].append([k, r_L, range_])
    d_a = {}
    for key, L in d_L.items():
        L.sort(key=lambda x: x[0])
        d_a[key] = np.array(L)
    d_L_r_L = {}
    for key, a in d_a:
        i = a[:, 2].argmax()
        if i != 0:
            r_L_low = np.sqrt(a[i-1, 1] * a[i, 1])
        else:
            r_L_low = -1
        if i != np.size(a) - 1:
            r_L_high = np.sqrt(a[i, 1] * a[i + 1, 1])
        else:
            r_L_high = -1
        d_L_r_L[key] = [r_L_low, r_L_high]
    return d_L_r_L
