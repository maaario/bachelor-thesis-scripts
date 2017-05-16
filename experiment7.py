# For one fi, plot approximated distribution with gaussians

import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import norm, kstest,shapiro

from experiment3 import generate_approximated_histograms
from kmerlight_theory import variance_of_absolute_estimate
from utils import *


mpl.rcParams.update({'font.size': 11})


def ecdf(x):
    xs = np.sort(x)
    ys = np.arange(1, len(xs)+1)/float(len(xs))
    return xs, ys


if __name__ == '__main__':
    # dirname_templ = 'experiment7/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}_fixed'  # 1000 t=1
    dirname_templ = 'experiment4/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}_fixed'  # 1000 t=7
    # dirname_templ = 'experiment7/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}_normal'  # 300 t=7
    # we must use expermental mean due to overestimation

    params = {'genome_size': 1000000, 'coverage': 50, 'error_rate': 0.02, 'k': 21, 'read_length': 100}
    # {'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.1, 'k': 21, 'read_length': 100}
    # {'genome_size': 1000000, 'coverage': 20, 'error_rate': 0.05, 'k': 21, 'read_length': 100}
    # {'genome_size': 1000000, 'coverage': 50, 'error_rate': 0.02, 'k': 21, 'read_length': 100}

    trials = 300
    t = 7
    r = 2**15

    dirname = dirname_templ.format(**params)
    generate_approximated_histograms(params, dirname, dirname, trials)

    exact_histogram = load_histograms_to_matrix(['{}/jellyfish.dist'.format(dirname)])[0]
    maxf = exact_histogram.shape[0] - 1
    approximated_histograms = load_histograms_to_matrix(
        ['{}/kmerlight_{:02d}.dist'.format(dirname, trial) for trial in range(trials)], maxf
    )

    all_abundances = sorted(range(1, 50), key=lambda x: exact_histogram[x])

    abundances = set(all_abundances[1::8])
    abundances = abundances.union({1, 2})
    abundances = sorted(abundances, key=lambda x: exact_histogram[x])

    h = int(np.sqrt(len(abundances)))
    w = int(np.ceil(len(abundances) / h))

    for ni, i in enumerate(abundances, start=1):
        fi_estimates = approximated_histograms[:, i]
        fi, F0 = exact_histogram[i], np.sum(exact_histogram)

        plt.figure('pdf')
        ax = plt.subplot(h, w, ni)
        ax.set_title('$f_{%d}=%d$' % (i, exact_histogram[i]))
        plt.yticks([])
        plt.tight_layout(h_pad=0.1, w_pad=0.1)
        plt.hist(fi_estimates, normed=True, bins=200)
        low, high = plt.xlim()
        x = np.arange(low, high, (high-low) / 1000)
        m, s = np.mean(fi_estimates), np.std(fi_estimates)
        plt.plot(x, norm.pdf(x, loc=m, scale=s), label='best normal approx.')
        m, s = np.mean(fi_estimates), np.sqrt(variance_of_absolute_estimate(F0, fi, r, t))
        plt.plot(x, norm.pdf(x, loc=m, scale=s), label='theoretical approx.')

        plt.figure('cdf')
        ax = plt.subplot(h, w, ni)
        ax.set_title('$f_{%d}=%d$' % (i, exact_histogram[i]))
        plt.tight_layout(h_pad=0.1, w_pad=0.1)
        plt.plot(*ecdf(fi_estimates))
        m, s = np.mean(fi_estimates), np.std(fi_estimates)
        plt.plot(x, norm.cdf(x,  loc=m, scale=s), label='best normal approx.')
        m, s = fi, np.sqrt(variance_of_absolute_estimate(F0, fi, r, t))
        plt.plot(x, norm.cdf(x,  loc=m, scale=s), label='theoretical approx.')

        #for trial in range(trials):
        #    m, s = fi, np.sqrt(variance_of_absolute_estimate(F0, fi_estimates[trial], r, t))
        #    plt.plot(x, norm.cdf(x, loc=m, scale=s))

    # KS test of normality
    # problem: distribution is discrete for i<20, i>35 ... variance is ok, distribution is discrete, not continuous
    plt.figure('KS test p-values')
    p_values_best = []
    p_values_est = []
    for i in all_abundances:
        fi_estimates = approximated_histograms[:, i]
        fi, F0 = exact_histogram[i], np.sum(exact_histogram)

        m, s = np.mean(fi_estimates), np.std(fi_estimates)
        _, P = kstest(fi_estimates, 'norm', args=(m, s), N=trials)
        # _, P = shapiro(fi_estimates)
        p_values_best.append(P)

        m, s = np.mean(fi_estimates), np.sqrt(variance_of_absolute_estimate(F0, fi, r, t))
        _, P = kstest(fi_estimates, 'norm', args=(m, s), N=trials)
        # _, P = shapiro(fi_estimates)
        p_values_est.append(P)
    plt.scatter(exact_histogram[all_abundances], p_values_best, label='best normal approx.', color='orange')
    plt.scatter(exact_histogram[all_abundances], p_values_est, label='theoretical approx.', color='green')
    plt.legend()
    plt.xscale('log')
    plt.hlines(0.05, *plt.xlim())

    plt.show()

