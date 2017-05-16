# Proof of concept - new kmerlight alternative based on maximizing the
# likelihood of observing t_i^(1), t_i^(2), ... t_i^(64)

import numpy as np


def new_fi_estimate(t_w, F0, r):
    w = np.arange(1, 20)
    p_cf = (1 - 1 / r) ** (F0 / 2 ** w - 1)
    c = 1 - np.sum(p_cf * 1/(2**w))  # probability that one k-mer dies
    # old version: np.sum(t_w) / -np.log(c)
    T = np.sum(t_w)
    return - T / (c - 1)

if __name__ == '__main__':
    F0 = 10 ** 7
    r = 2 ** 15

    w = np.arange(1, 20)
    p_cf = (1 / 2 ** w) * (1 - 1 / r) ** (F0 / 2 ** w - 1)

    fi = 100000
    t_w = np.round(p_cf * fi)
    print(t_w)
    print(new_fi_estimate(t_w, F0, r))