# calculate exact histogram
# calculate multiple approximate histograms
# plot mean, sd of histogram column relative errors in dependence from lambda
# difference from the article is that we will test it on different F_0

# is the fraction F_0 / f_i the only relevant factor in approximated histogram?
# no, the variance is a function of F0, f_i


import math
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from kmerlight_theory import variance_of_estimate_error
from utils import *

mpl.rcParams.update({'font.size': 13})


def generate_reads_and_exact_histogram(params, dirname):
    if not os.path.exists(dirname):
        os.makedirs(dirname)
        params['base_dir'] = dirname
        generate_reads(params)

        run(jellyfish_count.format(**params))
        run(jellyfish_hist.format(**params))

    else:
        print('Directory {} already exists, experiment data not generated anew.'.format(dirname), file=sys.stderr)


def generate_approximated_histograms(params, reads_dirname, dirname, trials=50, alternative=0):
    if not os.path.exists(reads_dirname):
        print('Directory {} with reads, does not exist.'.format(reads_dirname), file=sys.stderr)

    if not os.path.exists(dirname):
        os.makedirs(dirname)

    params['kmerlight_f'] = max(30, 5 * params['coverage'])
    params['reads_dir'] = reads_dirname
    params['base_dir'] = dirname
    params['kmerlight_alternative'] = alternative

    for trial in range(trials):
        filename = '{}/kmerlight_{:02d}.dist'.format(dirname, trial)
        if not os.path.exists(filename):
            run(run_kmerlight.format(**params))
            normalize_kmerlight_output('{}/kmerlight_partial.dist'.format(dirname), filename)


covered_fi = []
covered_F0 = []


def compute_approximation_error_statistics(dirname_exact, dirname_approximated, r, trials=50,):
    f = load_histogram('{}/jellyfish.dist'.format(dirname_exact))
    F0 = sum(f.values())
    for k, v in f.items():
        covered_F0.append(F0)
        covered_fi.append(v)

    hatf = []

    for trial in range(trials):
        hatf.append(load_histogram('{}/kmerlight_{:02d}.dist'.format(dirname_approximated, trial)))

    lambdas = []
    mean_errors = []
    sd_errors = []
    theory_sds = []

    for i in range(max(f.keys())):
        if f[i] != 0:
            lambdas.append(F0/f[i])

            errors = [(hatf[trial][i] - f[i])/f[i] for trial in range(trials)]
            mean_error = np.mean(errors)
            sd_error = np.std(errors)

            mean_errors.append(mean_error)
            sd_errors.append(sd_error)

            theory_sd = math.sqrt(variance_of_estimate_error(F0, f[i], r))
            theory_sds.append(theory_sd)

    s = sorted(zip(lambdas, mean_errors, sd_errors, theory_sds))
    lambdas = [s[i][0] for i in range(len(s))]
    mean_errors = [s[i][1] for i in range(len(s))]
    sd_errors = [s[i][2] for i in range(len(s))]
    theory_sds = [s[i][3] for i in range(len(s))]

    return lambdas, mean_errors, sd_errors, theory_sds


def plot_approximation_error_statistics(lambdas, mean_errors, sd_errors, theory_sds, label):
    plt.figure('mean of relative errors')
    plt.plot(lambdas, mean_errors, label='mean error ' + label)
    plt.xscale('log')
    plt.legend()

    plt.figure('sd. of relative errors')
    plt.plot(lambdas, sd_errors, label='error sd ' + label)
    plt.plot(lambdas, theory_sds, label='theory sd ' + label)
    plt.xscale('log')
    plt.legend()


def plot_approximation_error_statistics2(lambdas, mean_errors, sd_errors, theory_sds):
    plt.figure('statistics of relative errors')
    plt.scatter(lambdas, mean_errors, label='mean of relative errors')
    plt.scatter(lambdas, sd_errors, label='sd. of relative errors')
    plt.xscale('log')
    # plt.ylim(0.1, 8)
    # plt.yscale('log')
    plt.legend()
    plt.ylabel(r'statistic of $(\hat f_i - f_i) / f_i$')
    plt.xlabel(r'$\lambda = F_0 / f_i$')


def plot_sd_theory(lambdas1, sd_errors1, lambdas2, sd_errors2, theory_sds):
    lambdas1 = [1.3 * 10**7 / l for l in lambdas1]
    lambdas2 = [1.3 * 10 ** 7 / l for l in lambdas2]

    plt.scatter(lambdas1, sd_errors1, label='original Kmerlight')
    plt.scatter(lambdas2, sd_errors2, label='modified Kmerlight')
    plt.plot(lambdas1, theory_sds, label='theory sd', color='green')
    plt.xscale('log')
    # plt.ylim(0.1, 8)
    # plt.yscale('log')
    plt.legend()
    plt.ylabel(r'sd. of $(\hat f_i - f_i) / f_i$')
    plt.xlabel(r'$f_i$')
    # plt.xticks(1.18 * 10**np.arange(0, 6), [r'$10^{}$'.format(i) for i in range(7, 1, -1)])


def plot_hat_fi_distribution(dirname_exact, dirname_approximated, trials, i):
    f = load_histogram('{}/jellyfish.dist'.format(dirname_exact))
    F0 = sum(f.values())

    fi = f.get(i, 0)
    fi_s = []

    for trial in range(trials):
        hist = load_histogram('{}/kmerlight_{:02d}.dist'.format(dirname_approximated, trial))
        fi_s.append(hist.get(i, 0))

    plt.figure('hat_fi_distribution i={}'.format(i))
    plt.hist(fi_s)
    plt.vlines([fi], 0, plt.ylim()[1], color='red')


def plot_many_datasets_comparison():
    dirname_templ = 'experiment3/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}'
    exp2_dirname_templ = 'experiment2/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}/trial00'
    param_sets = [
        ({'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100}, dirname_templ),
        ({'genome_size': 100000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100}, dirname_templ),
        ({'genome_size': 1000000, 'coverage': 50, 'error_rate': 0.01, 'k': 21, 'read_length': 100}, dirname_templ),
        ({'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.1, 'k': 21, 'read_length': 100}, dirname_templ),
        ({'genome_size': 1000000, 'coverage': 0.5, 'error_rate': 0.01, 'k': 21, 'read_length': 100}, exp2_dirname_templ),
        ({'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.05, 'k': 21, 'read_length': 100}, exp2_dirname_templ),
    ]

    for params, reads_dirname_templ in param_sets:
        dirname = dirname_templ.format(**params)
        reads_dirname = reads_dirname_templ.format(**params)

        if params['coverage'] == 50:
            trials = 100
        else:
            trials = 50

        # generate_approximated_histograms(params, reads_dirname, dirname, trials)
        lambdas, mean_errors, sd_errors, theory_sds = compute_approximation_error_statistics(
            reads_dirname, dirname, 2**19, trials)
        plot_approximation_error_statistics(lambdas, mean_errors, sd_errors, theory_sds,
                                            'L{genome_size}_c{coverage}_e{error_rate}'.format(**params))

    plt.figure('covered fi and F0 values')
    plt.scatter(covered_fi, covered_F0)
    plt.xlabel('fi (histogram columns)')
    plt.ylabel('F0 (sum of columns)')
    plt.xscale('log')


def plot_one_dataset_comparison():
    dirname_templ = 'experiment4/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}_normal'
    params = {'genome_size': 1000000, 'coverage': 50, 'error_rate': 0.02, 'k': 21, 'read_length': 100}
    trials = 300
    dirname = dirname_templ.format(**params)
    lambdas, mean_errors, sd_errors, theory_sds = compute_approximation_error_statistics(
        dirname, dirname, 2**15, trials)
    # plot_approximation_error_statistics2(lambdas, mean_errors, sd_errors, theory_sds)

    dirname_templ = 'experiment4/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}_fixed'
    dirname = dirname_templ.format(**params)
    lambdas2, _, sd_errors2, _ = compute_approximation_error_statistics(
        dirname, dirname, 2 ** 15, trials)

    plt.figure('statistics of relative errors')
    plt.subplot(121)
    plot_sd_theory(lambdas, sd_errors, lambdas2, sd_errors2, theory_sds)
    plt.subplot(122)
    plot_sd_theory(lambdas, sd_errors, lambdas2, sd_errors2, theory_sds)


if __name__ == '__main__':
    #plot_many_datasets_comparison()
    plot_one_dataset_comparison()
    plt.show()

