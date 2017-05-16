# Compare kmerlight approximated histograms with fixed w* to normal kmerlight run
# Also compare to new kmerlight alternative - likelihood
# The data was generated using script experiment4.py

import matplotlib.pyplot as plt
import matplotlib as mpl

from kmerlight_theory import variance_of_absolute_estimate
from experiment4 import plot_covest_errors_boxplot
from utils import *

mpl.rcParams.update({'font.size': 12})

kmerlight_variants = ['normal', 'fixed']  # , 'like']
labels = ['original Kmerlight', 'modified Kmerlight']

dirname_templ = 'experiment4/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}_{variant}'

params = {'genome_size': 1000000, 'coverage': 50, 'error_rate': 0.02, 'k': 21, 'read_length': 100}      # 300 trials
# params = {'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.1, 'k': 21, 'read_length': 100}      # 50 trials
# params = {'genome_size': 1000000, 'coverage': 20, 'error_rate': 0.05, 'k': 21, 'read_length': 100}    # 50 trials

trials = 300
estimates = []

dirname = dirname_templ.format(variant='normal', **params)
f = load_histograms_to_matrix(['{}/jellyfish.dist'.format(dirname)])[0]
maxf = f.shape[0] - 1
F0 = np.sum(f)

for variant, label in zip(kmerlight_variants, labels):
    dirname = dirname_templ.format(variant=variant, **params)

    approximated_histograms = load_histograms_to_matrix(
        ['{}/kmerlight_{:02d}.dist'.format(dirname, trial) for trial in range(trials)], maxf
    )
    histogram_errors = approximated_histograms - f
    errors, sds = np.mean(histogram_errors, axis=0), np.std(approximated_histograms, axis=0)

    plt.figure('Means')
    plt.plot(errors, marker='o', label=label)
    plt.hlines(0, *plt.xlim())
    plt.xlabel('i')
    plt.ylabel(r'mean of $(\hat f_i - f_i)$')
    plt.legend()

    plt.figure('St. deviations')
    plt.plot(sds, marker='o', label=label)
    plt.xlabel('i')
    plt.ylabel(r'st. deviation of $(\hat f_i - f_i)$')

    # data were generated only fo 50 trials
    # estimate = load_coverage_estimates(
    #     ['{}/kmerlight_{:02d}.est'.format(dirname, trial) for trial in range(trials)])
    # estimates.append(estimate)

plt.figure('St. deviations')
sds_theory = [np.sqrt(variance_of_absolute_estimate(F0, f[i], 2 ** 15)) for i in range(len(f))]
plt.plot(sds_theory, label='theory')
plt.legend()

# plot_covest_errors_boxplot(np.vstack(estimates) - params['coverage'], labels)

plt.show()
