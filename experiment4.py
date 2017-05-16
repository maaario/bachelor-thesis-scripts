# At first used to compare covest results with and without correction (mean removal)
# dataset (c, e, L) = one genome, one reads file
#   generate many kmerlight histograms (as in exp 3)
#   calculate error statistics (mu, sigma for different i) for each kmerlight hist / average
#   calculate covest estimates (average/boxplot) without and with correction

# Currently used to plot best and worst approximated histograms for covest

import matplotlib.pyplot as plt

from experiment3 import generate_approximated_histograms
from utils import *


def run_covest_on_histograms(dirname, filename_templ, params, trials):
    mparams = dict(**params)
    mparams['histogram'] = 'jellyfish.dist'
    mparams['likelihood'] = 'normal'
    mparams['likelihood-data-file'] = 'experiments/error_stats.errs'

    mparams['base_dir'] = dirname
    run(run_covest.format(**mparams),
        output='{}/jellyfish.est'.format(dirname),
        errors='{}/jellyfish.est.log'.format(dirname))

    for trial in range(trials):
        filename = filename_templ.format(trial)
        filepath = dirname + '/' + filename
        mparams['histogram'] = filename + '.dist'
        run(run_covest.format(**mparams), output=filepath + '.est', errors=filepath + '.est.log')


def plot_histograms(histograms, labels, error_bars=None, cut_start=1, cut_end=None):
    n, maxf = histograms.shape
    if cut_end is None:
        cut_end = maxf

    gap = 0.2
    bar_width = (1 - gap) / n
    x = np.arange(cut_start, cut_end)

    for i in range(n):
        offset = - (1 - gap) / 2 + (i + 0.5) * bar_width
        plt.bar(x + offset, histograms[i, cut_start:cut_end], bar_width, label=labels[i])
        if error_bars is not None and error_bars[i] is not None:
            plt.vlines(x + offset,
                       histograms[i, cut_start:cut_end] - error_bars[i][cut_start:cut_end],
                       histograms[i, cut_start:cut_end] + error_bars[i][cut_start:cut_end])

    plt.xlabel("i")
    plt.xticks(np.arange(cut_start, cut_end,2))
    plt.ylabel(r"$f_i$", rotation='horizontal')
    plt.yscale('log')
    plt.legend()


def plot_covest_errors_boxplot(errors, labels, figname='Covest errors on kmerlight histogram'):
    colors = plt.rcParams['axes.color_cycle'][:len(errors)]

    plt.figure(figname)
    bp = plt.boxplot([x for x in errors], notch=True, widths=0.8, patch_artist=True, labels=labels)
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    plt.setp(bp['medians'], color='black')
    plt.axhline(y=0)


def remove_means(dirname, mean_errors, trials=50):
    for trial in range(trials):
        hist = load_histogram('{}/kmerlight_{:02d}.dist'.format(dirname, trial))
        for i in hist:
            hist[i] = max(0, hist[i] - int(mean_errors[i-1]))

        with open('{}/kmerlight_{:02d}_rm.dist'.format(dirname, trial), 'w') as f:
            for i, cnt in sorted(hist.items()):
                f.write('{} {}\n'.format(i, cnt))


if __name__ == '__main__':
    dirname_templ = 'experiment5.1/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}_normal'

    param_sets = [
        # {'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.1, 'k': 21, 'read_length': 100},
        # worst is t=42, f5 is 10^3 instead of 10^2 ... is this the problem?
        # yes, partly - if f5 is decreased to 10^2, estimate goes from 14 to 7

        # {'genome_size': 1000000, 'coverage': 20, 'error_rate': 0.05, 'k': 21, 'read_length': 100},
        # here enhancing the tail does not help, but reducing the overestimate at lower fi (4-5) does help
        # the effect of the mean is prevalent

        {'genome_size': 1000000, 'coverage': 50, 'error_rate': 0.02, 'k': 21, 'read_length': 100},
        # local changes do not help with worst - many columns are overestimated
    ]

    trials = 20

    for params in param_sets:
        dirname = dirname_templ.format(**params)
        # generate_reads_and_exact_histogram(params, dirname)
        # generate_approximated_histograms(params, dirname, dirname, trials=trials)
        # run_covest_on_histograms(dirname, 'kmerlight_{:02d}', params, trials=trials)

        exact_histogram = load_histograms_to_matrix(['{}/jellyfish.dist'.format(dirname)])
        maxf = exact_histogram.shape[1] - 1
        approximated_histograms = load_histograms_to_matrix(
            ['{}/kmerlight_{:02d}.dist'.format(dirname, trial) for trial in range(trials)], maxf
        )
        histogram_errors = approximated_histograms - exact_histogram
        mean_errors, sds = np.mean(histogram_errors, axis=0), np.std(approximated_histograms, axis=0)

        save_error_statistics(dirname, mean_errors, sds)

        # remove_means(dirname, mean_errors)
        # run_covest_on_histograms(dirname, 'kmerlight_{:02d}_rm', params, trials=50)

        coverage_estimates = load_coverage_estimates(
            ['{}/kmerlight_{:02d}.est'.format(dirname, trial) for trial in range(trials)])
        # coverage_estimates_rm = load_coverage_estimates(
        #     ['{}/kmerlight_{:02d}_rm.est'.format(dirname, trial) for trial in range(trials)])
        # errors = np.vstack([coverage_estimates, coverage_estimates_rm]) - params['coverage']
        # plot_covest_errors_boxplot(errors, [r'$\hat c_k - c$', 'mean error removed'])

        plt.figure('Exact vs. average approximated histogram')
        plot_histograms(
            np.vstack([exact_histogram, exact_histogram + mean_errors]), ['jellyfish', '50 x kmerlight'], [None, sds]
        )

        abs_estimate_errors = np.abs(coverage_estimates - params['coverage'])
        best = np.argmin(abs_estimate_errors)
        worst = np.argmax(abs_estimate_errors)
        good_dist = np.percentile(abs_estimate_errors, 20)
        bad_dist = np.percentile(abs_estimate_errors, 80)

        plt.figure('Covest estimates on kmerlight histograms')
        plt.scatter(range(trials), coverage_estimates)
        plt.xlabel('trial')
        plt.scatter(best, coverage_estimates[best], color='green')
        plt.scatter(worst, coverage_estimates[worst], color='red')
        c = params['coverage']
        lims = plt.xlim()
        plt.hlines(c, *lims)
        plt.hlines([c - good_dist, c + good_dist], *lims, color='green')
        plt.hlines([c - bad_dist, c + bad_dist], *lims, color='red')

        plt.figure('Exact vs. worst vs. best histogram')
        plot_histograms(
            np.vstack([exact_histogram, approximated_histograms[worst], approximated_histograms[best]]),
            labels=['exact', 'worst (t={})'.format(worst), 'best (t={})'.format(best)]
        )

        plt.figure('Exact vs. bad vs. good histograms')
        good = approximated_histograms[abs_estimate_errors <= good_dist, :]
        bad = approximated_histograms[abs_estimate_errors >= bad_dist, :]
        plot_histograms(
            np.vstack([exact_histogram, np.mean(good, axis=0), np.mean(bad, axis=0)]),
            labels=['exact', 'bad', 'good'],
            error_bars=[None, np.std(bad, axis=0), np.std(good, axis=0)]
        )

        plt.show()
