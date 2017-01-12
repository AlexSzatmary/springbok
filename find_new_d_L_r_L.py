from collections import defaultdict
import glob
import numpy as np


def argmax_last(a):
    return np.size(a) - 1 - np.argmax(a[::-1])


def argmax_end_of_first_peak(b):
    if np.sum(b == np.max(b)) > 1:
        # Have to do it the hard way
        i = np.argmax(b)
        while True:
            if b[i + 1] == b[i] and i < np.size(b):
                i += 1
            else:
                break
    else:
        i = np.argmax(b)
    return i


def find(a):
    d_L = defaultdict(lambda: [])
    for row in a:
        (gamma_L, phi_E, k, r_L, range_) = row
        d_L[gamma_L, phi_E].append([k, r_L, range_])
    d_a = {}
    for key, L in d_L.items():
        L.sort(key=lambda x: x[1])
        d_a[key] = np.array(L)
    d_L_r_L = {}
    for key, a in d_a.items():
        i = argmax_end_of_first_peak(a[:, 2])
        if i > 1:
            r_L_low = np.sqrt(a[i-1, 1] * a[i, 1])
        elif i == 1:
            r_L_low = a[i, 1] / 10
        else:
#            r_L_low = -1
            r_L_low = a[i + 1, 1] / 100
        if i != np.size(a, 0) - 1:
#            print([key, a, i])
            if i != 0:
                r_L_high = np.sqrt(a[i, 1] * a[i + 1, 1])
            else:
                r_L_high = a[i + 1, 1] / 10
        else:
            r_L_high =a[i, 1] * 10
        d_L_r_L[key] = [r_L_low, r_L_high]
    return (d_a, d_L_r_L)
