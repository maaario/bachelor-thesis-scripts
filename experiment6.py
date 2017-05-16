# Test functionality of different error functions for covest

# normal histogram
#   -> apply specific modification (generate "estimate" with different variances for different i) many times
#       -> use covest + covest with different error functions

# generate modified_xx.dist, jellyfish.est, modified_xx_covesttype.est from jellyfish.dist

import matplotlib.pyplot as plt

from experiment4 import plot_covest_errors_boxplot
from utils import *


def save_stds(filename, stds):
    with open(filename, 'w') as f:
        f.write(' '.join(map(str, stds)))
        f.write('\n')


def modify_histogram(histogram, stds):
    return np.random.normal(histogram, stds).clip(min=0)


def save_histogram_from_array(filename, histogram):
    with open(filename, 'w') as f:
        f.write('\n'.join([' '.join(map(str, map(int, x))) for x in enumerate(histogram[1:], start=1)]))


def generate_modified_histograms(dirname, exact_histogram, stds, trials=50):
    for trial in range(trials):
        filename = '{}/modified_{:02d}.dist'.format(dirname, trial)
        if not os.path.exists(filename):
            histogram = modify_histogram(exact_histogram, stds)
            save_histogram_from_array(filename, histogram)


def run_covest_on_histograms(dirname, params, trials, covest_type='normal'):
    mparams = dict(**params)
    mparams['base_dir'] = dirname
    mparams['likelihood'] = covest_type
    mparams['likelihood-data-file'] = '{}/stds.txt'.format(dirname)

    mparams['histogram'] = 'jellyfish.dist'
    filename = '{}/jellyfish_{}.est'.format(dirname, covest_type)
    if not os.path.exists(filename):
        run(run_covest_alternative.format(**mparams), output=filename, errors=filename + '.log')

    for trial in range(trials):
        histogram_filename = 'modified_{:02d}.dist'.format(trial)
        mparams['histogram'] = histogram_filename
        estimate_filename = '{}/modified_{:02d}_{}.est'.format(dirname, trial, covest_type)
        if not os.path.exists(estimate_filename):
            run(run_covest_alternative.format(**mparams), output=estimate_filename, errors=estimate_filename + '.log')

if __name__ == '__main__':
    trials = 100
    params = {'genome_size': 1000000, 'coverage': 50, 'error_rate': 0.02, 'k': 21, 'read_length': 100}

    dirname = 'experiment6/normal_sf/'
    # dirname = 'experiment6/normal_sconst10^5/'
    # dirname = 'experiment6/normal_ssqrtf/'
    # dirname = 'experiment6/normal_sf2/'
    # dirname ='experiment6/normal_srandom2f'
    # dirname ='experiment6/normal_s0.33f'


    exact_histogram = load_histograms_to_matrix(['{}/jellyfish.dist'.format(dirname)])[0]
    # stds = exact_histogram / 3  # /normal_s0.33f
    # stds = [np.random.random() * fi for fi in exact_histogram]
    # stds = [fi for fi in exact_histogram]

    # save_stds('{}/stds.txt'.format(dirname), stds)
    # generate_modified_histograms(dirname, exact_histogram, stds, trials)
    modified_histograms = load_histograms_to_matrix(
        ['{}/modified_{:02d}.dist'.format(dirname, trial) for trial in range(trials)])

    methods = ['normal',
               # 'minsq',
               'wminsq']

    errors = []
    for method in methods:
        # run_covest_on_histograms(dirname, params, trials, covest_type=method)
        coverage_estimates = load_coverage_estimates(
            ['{}/modified_{:02d}_{}.est'.format(dirname, trial, method) for trial in range(trials)])

        plt.figure('Covest estimates on modified histograms per trial')
        plt.scatter(range(trials), coverage_estimates - params['coverage'], label=method)

        errors.append(coverage_estimates - params['coverage'])

    plt.figure('Covest estimates on modified histograms per trial')
    plt.ylim(-params['coverage'], params['coverage'])
    plt.legend()

    plt.figure('Covest estimates on modified histograms boxplot')
    plot_covest_errors_boxplot(errors, ['log-likelihood', 'w. min. sq. err.'], figname='Covest estimates on modified histograms boxplot')
    plt.ylim(-params['coverage'], params['coverage'])

    plt.show()
