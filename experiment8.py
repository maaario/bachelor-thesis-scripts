# Compare log. likelihood and sum of square errors

import matplotlib.pyplot as plt

from experiment6 import modify_histogram, save_histogram_from_array, save_stds
from utils import *

params = {'genome_size': 1000000, 'k': 21, 'read_length': 100, 'likelihood-data-file': 'error_stats.errs'}


def compare_covest_types(histograms_dir, output_dir, covest_types, default_params, generate=False):
    params = dict(**default_params)
    params['base_dir'] = histograms_dir
    est_files = [[] for _ in covest_types]

    for coverage in [0.1, 0.5, 1, 2, 4, 10, 50]:
        for error_rate in [0.01, 0.03, 0.05, 0.1]:
            params['coverage'] = coverage
            params['error_rate'] = error_rate

            for trial in range(5):
                filename = 'simulated{trial}_c{coverage}_e{error_rate}_k{k}'.format(**params, trial=trial)
                params['histogram'] = filename + '.dist'

                for i, ct in enumerate(covest_types):
                    est_filename = '{}/{}_{}.est'.format(output_dir, filename, ct)
                    est_files[i].append(est_filename)
                    if generate:
                        params['likelihood'] = ct
                        if ct == 'wminsq':
                            params['likelihood-data-file'] = '{}/{}.dist.errs'.format(histograms_dir, filename)
                        run(run_covest_alternative.format(**params), output=est_filename, errors=est_filename + '.log')

    means = []
    stds = []

    for i, ct in enumerate(covest_types):
        estimates = load_coverage_estimates(est_files[i])
        estimates = estimates.reshape([len(estimates) // 5, 5])
        means.append(np.mean(estimates, axis=1))
        stds.append(np.std(estimates, axis=1))

    return means, stds

all_coverages = np.repeat([0.1, 0.5, 1, 2, 4, 10, 50], 4)

def plot_comparison(means, stds, covest_types):
    datapoints = len(means[0])
    x = list(range(datapoints))

    plt.figure('Means')
    for i, ct in enumerate(covest_types):
        plt.scatter(x, (means[i] - all_coverages), label=ct)
    # plt.yscale('log')
    # plt.ylim(0,5)
    plt.vlines(np.arange(0, datapoints+1, 4) - 0.5, *plt.ylim())
    plt.hlines(0, *plt.xlim())
    # plt.hlines([0.1, 0.5, 1, 2, 4, 10, 50], np.arange(0, datapoints+1, 4) - 0.5, np.arange(0, datapoints+1, 4) + 3.5)
    plt.xticks([])
    plt.legend()

    plt.figure('Stdevs')
    for i, ct in enumerate(covest_types):
        plt.scatter(x, stds[i], label=ct)
    #plt.ylim(0, 7)
    plt.vlines(np.arange(0, datapoints+1, 4) - 0.5, *plt.ylim())
    plt.legend()

    plt.show()


def generate_modified_histograms(original_dirname, output_dirname, stdf):
    for name in os.listdir(original_dirname):
        exact_histogram = load_histograms_to_matrix(['{}/{}'.format(original_dirname, name)])[0]
        stds = np.array([stdf(fi) for fi in exact_histogram])
        histogram = modify_histogram(exact_histogram, stds)
        save_histogram_from_array('{}/{}'.format(output_dirname, name), histogram)
        save_stds('{}/{}.errs'.format(output_dirname, name), stds)

if __name__ == '__main__':
    # covest_types = ['normal', 'minsqlin']
    # means, stds = compare_covest_types(
    #     'experiment8/simulated_histograms', 'experiment8/compare_minsq_to_log',covest_types, params, generate=False
    # )
    # plot_comparison(means, stds, covest_types)

    # # generate_modified_histograms(
    # #    'experiment8/simulated_histograms', 'experiment8/simulated_histograms_s0.3f', lambda x: 0.3 * x)
    # covest_types = ['normal', 'wminsq']
    # means, stds = compare_covest_types(
    #     'experiment8/simulated_histograms_s0.3f', 'experiment8/compare_wminsq_to_log_s0.3f',
    #     covest_types, params, generate=False
    # )
    # plot_comparison(means, stds, covest_types)

    # generate_modified_histograms(
    #     'experiment8/simulated_histograms', 'experiment8/simulated_histograms_srandomf',
    #     lambda x: np.random.random() * x)
    covest_types = ['normal', 'wminsq']
    means, stds = compare_covest_types(
        'experiment8/simulated_histograms_srandomf', 'experiment8/compare_wminsq_to_log_srandomf',
        covest_types, params, generate=False
    )
    plot_comparison(means, stds, covest_types)
