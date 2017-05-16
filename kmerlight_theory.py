import numpy as np
import matplotlib.pyplot as plt


def calc_w_star(F0, r):
    # return max(1, np.ceil(np.log(- np.log(1 - 1 / r) * F0) / np.log(2)))  # result of differentiation
    return np.argmax([(1 - 1 / r) ** (F0 / (2 ** w) - 1) / 2 ** w for w in range(1, 60)]) + 1


def calc_p_sample(F0, r):
    w_star = calc_w_star(F0, r)

    # collisions from all k-mers, little bit more precise, worse for calculations
    # return 1 / (2 ** w_star) * (1 - 1 / (2**w_star * r)) ** (F0 - 1)

    return 1 / (2 ** w_star) * (1 - 1 / r) ** (F0 / 2 ** w_star - 1)


def variance_of_estimate_error(F0, fi, r, t=7):
    # Var( (hat fi - fi)/fi ),  hat fi ~ 1/p * Bin(fi, P(kmer survives at level w*))
    p_sample = calc_p_sample(F0, r)
    medians_factor = 1 / (4 * t / (2 * np.pi))
    variance = medians_factor * (1 - p_sample) / (p_sample * fi)
    return variance


def variance_of_absolute_estimate(F0, fi, r, t=7):
    # Var(hat fi),  hat fi ~ 1/p * Bin(fi, P(kmer survives at level w*))
    p_sample = calc_p_sample(F0, r)
    if t == 1:
        medians_factor = 1
    else:
        medians_factor = 1 / (4 * t / (2 * np.pi))
    variance = medians_factor * fi * (1 - p_sample) / p_sample
    return variance

if __name__ == '__main__':
    F0 = 1.2 * 10 ** 7
    r = 2 ** 15

    plt.figure('Discrete and continuous choice of best w*')
    w = np.arange(1, 20)
    p = (1 - 1 / r) ** (F0 / 2 ** w - 1)
    plt.scatter(w, p, label='this should be in (1/4, 1/2)')

    start = np.log(- np.log(1 - 1 / r) * F0) / np.log(2)
    w = np.arange(start, 20, 0.1)
    dp = 1 + np.log(1 - 1 / r) * F0 / 2**w
    plt.plot(w, dp, label='derivative of $p_{cf}$ = 0')

    plt.hlines([0, 1 / 2, 1 / 4], *plt.xlim())
    plt.legend()

    import matplotlib as mpl
    mpl.rcParams.update({'font.size': 12})

    plt.figure('Empty, collision-free, collision counters')
    w = np.arange(1, 17)
    p_empty = (1 - 1 / r) ** (F0 / 2**w)
    p_cf = (F0 / 2**w) * (1 / r) * (1 - 1 / r) ** (F0 / 2**w - 1)
    p_c = 1 - p_empty - p_cf

    w_cont = np.arange(1, 17, 0.1)
    p_empty_cont = (1 - 1 / r) ** (F0 / 2 ** w_cont)
    p_cf_cont = (F0 / 2 ** w_cont) * (1 / r) * (1 - 1 / r) ** (F0 / 2 ** w_cont - 1)
    p_c_cont = 1 - p_empty_cont - p_cf_cont

    plt.scatter(w, p_empty, label='empty', color='blue')
    plt.plot(w_cont, p_empty_cont, color='blue')
    plt.scatter(w, p_cf, label='collision-free', color='green')
    plt.plot(w_cont, p_cf_cont, color='green')
    plt.scatter(w, p_c, label='collisions', color='red')
    plt.plot(w_cont, p_c_cont, color='red')

    plt.xticks(w)
    plt.xlabel('w')
    plt.ylabel('fractions of counters at one level')
    plt.legend()

    plt.figure('Probability of survival of 1 k-mer')
    fi = 2 ** np.arange(20)

    p_one_dies = np.prod((p_c / (p_cf + p_c)) ** (1 / 2**w))
    p_at_least_one_survives = 1 - p_one_dies**fi
    plt.scatter(fi, p_at_least_one_survives)

    p_one_dies = 1 - np.sum((1 / 2**w) * (1 - 1 / (2**w * r)) ** (F0 - 1) )
    p_at_least_one_survives = 1 - p_one_dies**fi
    plt.scatter(fi, p_at_least_one_survives)

    plt.xscale('log')
    plt.show()
