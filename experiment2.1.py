# As experiment 2 but this time compare covest, kmerlight and modified kmerlight

import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import kstest, ttest_rel


from experiment2 import generate_dataset, get_coverage_estimates, group_boxplot
from kmerlight_theory import calc_w_star
from utils import *

mpl.rcParams.update({'font.size': 13})


dirname_temp = 'experiment2.1/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}'
trial_dirname_temp = dirname_temp + '/trial{trial:02d}'


def generate_table(data, titles):
    # data[dataset][algorithm used][trial] = estimate

    for i, dataset in enumerate(data):
        means = []
        stds = []

        for trials in dataset:
            _trials = [t for t in trials if t != -50]
            means.append(np.mean(_trials))
            stds.append(np.std(_trials))

        print(('{}' + ' & {:.3f}'*6 + '\\\\ \hline').format(titles[i], *(means+stds)))


def differences_caused_by_approxiamtion_vs_param_sets(param_sets, titles, trials):
    covest_estimates = []
    kmerlight_coverage_estimates = []
    mkmerlight_coverage_estimates = []

    coverage_estimate_differences = []
    coverage_estimate_differences2 = []
    coverage_estimate_differences3 = []

    for params in param_sets:
        generate_dataset(params, trial_dirname_temp, trials)

        # first_trial_dir = trial_dirname_temp.format(**params, trial=0)
        # exact_histogram = load_histograms_to_matrix(['{}/jellyfish.dist'.format(first_trial_dir)])[0]
        # F0 = np.sum(exact_histogram)
        # r = 2**15
        # wplus = calc_w_star(F0, r)
        # w = np.array([wplus-1, wplus, wplus+1])
        # p_cf = (F0 / 2 ** w) * (1 / r) * (1 - 1 / r) ** (F0 / 2 ** w - 1)
        # print(p_cf)

        jellyfish_coverage, kmerlight_coverage, mkmerlight_coverage = np.array(get_coverage_estimates(
            dirname_temp.format(**params), ['jellyfish', 'kmerlight', 'mkmerlight']
        ))

        # for dataset in [jellyfish_coverage, kmerlight_coverage, mkmerlight_coverage]:
        #     _, p = kstest(dataset, 'norm', args=(np.mean(dataset), np.std(dataset)))
        #     print(p)

        x = kmerlight_coverage - mkmerlight_coverage
        _, KSp = kstest(x, 'norm', args=(np.mean(x), np.std(x)))
        _, STp = ttest_rel(kmerlight_coverage, mkmerlight_coverage)
        print(KSp, STp)

        covest_estimates.append(jellyfish_coverage - params['coverage'])
        kmerlight_coverage_estimates.append(kmerlight_coverage - params['coverage'])
        mkmerlight_coverage_estimates.append(mkmerlight_coverage - params['coverage'])

        # if covest has larger variance than variance introduced by kmerlight, than we must consider differences
        # since estimates are corellated
        coverage_estimate_differences.append(kmerlight_coverage - jellyfish_coverage)
        coverage_estimate_differences2.append(mkmerlight_coverage - jellyfish_coverage)
        coverage_estimate_differences3.append(kmerlight_coverage - mkmerlight_coverage)

    data = list(zip(
        covest_estimates,
        kmerlight_coverage_estimates,
        mkmerlight_coverage_estimates,
        # coverage_estimate_differences,
        # coverage_estimate_differences2,
        # coverage_estimate_differences3
    ))

    labels = [
        r'$\hat c_j - c$',
        r'$\hat c_{ok} - c$',
        r'$\hat c_{mk} - c$',
        # r'$\hat c_{ok} - \hat c_j$',
        # r'$\hat c_{mk} - \hat c_j$',
        # r'$\hat c_{ok} - \hat c_{mk}$',
    ]

    generate_table(data, titles)
    group_boxplot(data, labels, titles)
    # plt.tight_layout()
    plt.show()
    plt.clf()

if __name__ == '__main__':
    trials = 50

    titles = ['e=0', 'e=0.01', 'e=0.05', 'e=0.1']
    param_sets = [
        {'genome_size': 1000000, 'coverage': 10, 'error_rate': 0, 'k': 21, 'read_length': 100},
        {'genome_size': 1000000, 'coverage': 10, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        {'genome_size': 1000000, 'coverage': 10, 'error_rate': 0.05, 'k': 21, 'read_length': 100},
        {'genome_size': 1000000, 'coverage': 10, 'error_rate': 0.1, 'k': 21, 'read_length': 100},
    ]
    differences_caused_by_approxiamtion_vs_param_sets(param_sets, titles, trials)

    titles = ['c=0.5', 'c=2', 'c=10', 'c=50']
    param_sets = [
        {'genome_size': 1000000, 'coverage': 0.5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        {'genome_size': 1000000, 'coverage': 2,   'error_rate': 0.01, 'k': 21, 'read_length': 100},
        # {'genome_size': 1000000, 'coverage': 5,   'error_rate': 0.01, 'k': 21, 'read_length': 100},
        {'genome_size': 1000000, 'coverage': 10,  'error_rate': 0.01, 'k': 21, 'read_length': 100},
        # {'genome_size': 1000000, 'coverage': 20, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        {'genome_size': 1000000, 'coverage': 50,  'error_rate': 0.01, 'k': 21, 'read_length': 100},
    ]
    differences_caused_by_approxiamtion_vs_param_sets(param_sets, titles, trials)

    titles = [r'$L=10^5$', r'$L=10^6$', r'$L=10^7$']
    param_sets = [
        {'genome_size': 100000, 'coverage': 10, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        {'genome_size': 1000000, 'coverage': 10, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        {'genome_size': 10000000, 'coverage': 10, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
    ]
    differences_caused_by_approxiamtion_vs_param_sets(param_sets, titles, trials)
