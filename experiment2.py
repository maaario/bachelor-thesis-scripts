# which genome characteristics cause approximate histogram to influence covest?
# same comparison as experiment1, but each datapoint comes from 50 trials

# plot covest mean + variance (boxplot) with real/approximated histograms
# boxplots don't work due to correlation, use boxplots of differences


import sys

import matplotlib.pyplot as plt

from utils import *


def generate_dataset(params, dirname_tmpl, trials=50):
    for trial in range(trials):
        dirname = dirname_tmpl.format(**params, trial=trial)

        if not os.path.exists(dirname):
            os.makedirs(dirname)
            params['base_dir'] = dirname

            generate_reads(params)

            calculate_jellyfish_histogram(params)

            params['kmerlight_alternative'] = 0
            calculate_kmerlight_histogram(params)

            params['kmerlight_alternative'] = 1
            calculate_kmerlight_histogram(params, outfile='mkmerlight')

            run_covest_on_histograms(params, ['jellyfish', 'kmerlight', 'mkmerlight'])
        else:
            print('Directory {} already exists, experiment data not generated anew.'.format(dirname), file=sys.stderr)


def get_coverage_estimates(dirname, alternatives):
    coverages = [[] for _ in alternatives]

    for name in sorted(os.listdir(dirname)):
        subdir = os.path.join(dirname, name)
        if os.path.isdir(subdir) and os.path.exists('{}/{}.est'.format(subdir, alternatives[-1])):
            for i, alt in enumerate(alternatives):
                coverages[i].append(float(load_covest_results('{}/{}.est'.format(subdir, alt)).get('coverage', 0)))

    return coverages


def plot_boxplot(jellyfish_coverage, kmerlight_coverage):
    labels = ['jellyfish', 'kmerlight']
    plt.boxplot([jellyfish_coverage, kmerlight_coverage], labels=labels, notch=True, widths=0.8)
    plt.ylabel('coverage estimates')
    plt.title('Covest variance on exact and estimated histograms')
    plt.show()


def plot_scatter(jellyfish_coverage, kmerlight_coverage):
    plt.scatter(jellyfish_coverage, kmerlight_coverage)
    plt.axes().set_aspect('equal')

    lims = [
        min([*plt.xlim(), *plt.ylim()]),  # min of both axes
        max([*plt.xlim(), *plt.ylim()]),  # max of both axes
    ]
    plt.plot(lims, lims, 'k-')
    plt.xlim(*lims)
    plt.ylim(*lims)

    plt.xlabel('coverage est. on exact histogram')
    plt.ylabel('coverage est. on approximated histogram')
    plt.title('Covest estimates on exact and approximated histograms')
    plt.show()


def plot_difference(jellyfish_coverage, kmerlight_coverage):
    n = len(jellyfish_coverage)
    plt.scatter(list(range(n)), [kmerlight_coverage[i] - jellyfish_coverage[i] for i in range(n)])
    plt.axhline(y=0)
    plt.ylabel('coverage estimate differences')
    plt.xlabel('trial')
    plt.title('Covest estimate variance added by histogram approximation')
    plt.show()


def boxplot_one_param_set(data, labels, title):
    colors = plt.rcParams['axes.color_cycle'][:len(data)]
    bp = plt.boxplot(data, notch=True, widths=0.8, patch_artist=True, labels=labels)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    plt.setp(bp['medians'], color='black')
    plt.axhline(y=0)
    plt.title(title)


def group_boxplot(data, labels, titles):
    groups = len(data)
    for group in range(groups):
        plt.subplot(1, groups, group + 1)
        boxplot_one_param_set(data[group], labels, titles[group])


def differences_caused_by_approxiamtion_vs_param_sets(param_sets, titles, trials, figname):
    covest_estimates = []
    coverage_estimate_differences = []

    for params in param_sets:
        dirname_temp = 'experiment2/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}'
        trial_dirname_temp = dirname_temp + '/trial{trial:02d}'
        generate_dataset(params, trial_dirname_temp, trials)

        jellyfish_coverage, kmerlight_coverage = get_coverage_estimates(dirname_temp.format(**params),
                                                                        ['jellyfish', 'kmerlight'])
        # plot_boxplot(jellyfish_coverage, kmerlight_coverage)
        # data are corellated, boxplot is not useful
        # plot_scatter(jellyfish_coverage, kmerlight_coverage)
        # plot_difference(jellyfish_coverage, kmerlight_coverage)

        coverage = [jellyfish_coverage[i] - params['coverage'] for i in range(len(jellyfish_coverage))]
        covest_estimates.append(coverage)

        differences = [kmerlight_coverage[i] - jellyfish_coverage[i] for i in range(len(kmerlight_coverage))]
        coverage_estimate_differences.append(differences)

    group_boxplot(
        list(zip(covest_estimates, coverage_estimate_differences)),
        [r'$\hat c - c$', r'$\hat c_k - \hat c$'],
        titles,
    )
    plt.savefig('experiment2/' + figname)
    plt.show()
    plt.clf()


if __name__ == '__main__':
    trials = 50
    titles = ['default', 'c0.5', 'L100000']
    param_sets = [
        {'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        {'genome_size': 1000000, 'coverage': 0.5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        {'genome_size': 100000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
    ]
    differences_caused_by_approxiamtion_vs_param_sets(param_sets, titles, trials, 'unimportant_parameters.pdf')
    # length of genome and coverage are not that relevant

    titles = ['e0.01', 'e0.05', 'e0.1']
    param_sets = [
        {'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        {'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.05, 'k': 21, 'read_length': 100},
        {'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.1, 'k': 21, 'read_length': 100},
    ]
    differences_caused_by_approxiamtion_vs_param_sets(param_sets, titles, trials, 'error_rate.pdf')
    # covest overestimates coverage with higher error rate, with approximated histogram
    # probably because high error rate causes F0/f_i coefficient to be high -> inaccurate histograms
    # with higher error rate covest has also higher variance - due to variance in histogram
