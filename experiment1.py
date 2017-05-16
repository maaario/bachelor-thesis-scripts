# which genome characteristics cause approximate histogram to influence covest?
# generate various datasets
#   4 sets: different e, L, c, k
# compare histogram approximation errors
# plot graphs:
#   approximated histogram relative errors vs. parameter and column number
#   approximated vs. exact histogram (side by side + relative errors)
# compare covest results
#   covest results on approximated and exact histogram vs. parameter

import sys

import matplotlib.pyplot as plt

from utils import *


def generate_experimental_dataset(param_name, param_range, dirname_tmpl):
    params = dict(**default_params)
    
    for altered_param in param_range: 
        params[param_name] = altered_param
        dirname = dirname_tmpl.format(**params)
        params['base_dir'] = dirname

        if not os.path.exists(dirname):
            os.makedirs(dirname)
            generate_reads(params)
            calculate_jellyfish_histogram(params)
            calculate_kmerlight_histogram(params)
            run_covest_on_histograms(params, ['jellyfish', 'kmerlight'])
        else:
            print('Directory {} already exists, experiment data not generated anew.'.format(dirname), file=sys.stderr)


def plot_relative_histogram_approximation_errors_for_dataset(param_name, param_range, dirname_tmpl):
    errors = []
    params = dict(**default_params)

    for altered_param in param_range: 
        params[param_name] = altered_param
        dirname = dirname_tmpl.format(**params)
        params['base_dir'] = dirname

        if not os.path.exists(dirname):
            print('Directory {} not found, skip plotting relative errors.'.format(dirname), file=sys.stderr)
            errors.append([])
            continue

        jellyfish = load_histogram('{base_dir}/jellyfish.dist'.format(**params))
        kmerlight = load_histogram('{base_dir}/kmerlight.dist'.format(**params))

        all_kmers = sum(jellyfish.values())
        kmers_seen = 0

        max_i = max(jellyfish)
        for i in range(1, max_i+1):
            rel_error = (kmerlight[i] - jellyfish[i]) / jellyfish[i]
            errors.append((altered_param, rel_error, i))

            kmers_seen += jellyfish[i]
            if kmers_seen / all_kmers > 0.999:
                break

    # plot error distributions
    x = [x for x, y, z in errors]
    y = [y for x, y, z in errors]
    z = [z for x, y, z in errors]
    
    param_to_set = {param: i for i, param in enumerate(param_range)}
    equally_spaced_x = [param_to_set[l] for l in x] 
    
    plt.scatter(equally_spaced_x, y)
    for x, y, z in zip(equally_spaced_x, y, z):
        plt.annotate(z, xy=(x, y), xytext=(-5, -5), textcoords='offset points', ha='right', va='bottom')
    plt.xticks(range(len(param_range)), param_range)

    plt.xlabel(param_name)
    plt.ylabel("relative errors made by kmerlight")
    plt.title("Kmerlight relative errors vs. {}".format(param_name))
    plt.show()


def plot_coverage_estimates_for_dataset(param_name, param_range, dirname_tmpl, relative_error=True):
    jellyfish_coverage = []
    kmerlight_coverage = []

    params = dict(**default_params)

    for altered_param in param_range: 
        params[param_name] = altered_param
        dirname = dirname_tmpl.format(**params)
        params['base_dir'] = dirname
        
        if not os.path.exists(dirname):
            print('Directory {} not found, skip plotting coverage_estimate.'.format(dirname), file=sys.stderr)
            jellyfish_coverage.append(0)
            kmerlight_coverage.append(0)
            continue

        jellyfish_coverage.append(float(load_covest_results('{base_dir}/jellyfish.est'.format(**params)).get('coverage', 0)))
        kmerlight_coverage.append(float(load_covest_results('{base_dir}/kmerlight.est'.format(**params)).get('coverage', 0)))

    # plot coverage estimates
    equally_spaced_x = list(range(len(param_range)))
    
    if not relative_error:
        plt.scatter(equally_spaced_x, jellyfish_coverage, label='jellyfish')
        plt.scatter(equally_spaced_x, kmerlight_coverage, label='kmerlight')
        plt.xticks(equally_spaced_x, param_range)
        plt.xlabel(param_name)
        plt.ylabel("coverage estimated by covest")
        plt.title("Coverage estimates on kmerlight and jellyfish histograms")
        plt.legend()
        plt.show()
    else:
        rel_errors = []
        for i in range(len(param_range)):
            j = jellyfish_coverage[i]
            k = kmerlight_coverage[i]
            rel_errors.append((k-j)/j if j != 0 else 1)
        
        plt.scatter(equally_spaced_x, rel_errors)
        plt.xticks(equally_spaced_x, param_range)
        plt.xlabel(param_name)
        plt.ylabel("relative error of estimate")
        plt.title("Coverage estimate rel. errors on kmerlight vs. jellyfish histograms")
        plt.show()


def covest_estimate_on_kmerlight_and_jellyfish(param_name, param_range, dirname_tmpl):
    generate_experimental_dataset(param_name, param_range, dirname_tmpl)
    # plot_relative_histogram_approximation_errors_for_dataset(param_name, param_range, dirname_tmpl)

    plot_coverage_estimates_for_dataset(param_name, param_range, dirname_tmpl, relative_error=False)
    plot_coverage_estimates_for_dataset(param_name, param_range, dirname_tmpl, relative_error=True)


def plot_histograms(dirname):
    if not os.path.exists(dirname):
        print('Directory {} does not exist. Cannot display graph.'.format(dirname), file=sys.stderr)
        return

    jellyfish = load_histogram('{}/jellyfish.dist'.format(dirname))
    kmerlight = load_histogram('{}/kmerlight.dist'.format(dirname))

    plt.bar([x - 0.4 for x in jellyfish.keys()], jellyfish.values(), 0.4, label='jellyfish')
    plt.bar([x for x in kmerlight.keys()], kmerlight.values(), 0.4, label='kmerlight')
    plt.xlabel("i")
    plt.ylabel("H[i]")
    plt.yscale('log')
    plt.legend()
    plt.show()


def plot_histogram_approximation_errors(dirname, relative=True):
    if not os.path.exists(dirname):
        print('Directory {} does not exist. Cannot display graph.'.format(dirname), file=sys.stderr)
        return

    jellyfish = load_histogram('{}/jellyfish.dist'.format(dirname))
    kmerlight = load_histogram('{}/kmerlight.dist'.format(dirname))

    max_i = max(max(jellyfish.keys()), max(kmerlight.keys()))

    rel_errors = []
    for i in range(1, max_i+1):
        j = jellyfish[i]
        k = kmerlight[i]
        if relative:
            err = (k-j) / j if j != 0 else 1
        else:
            err = abs(k-j)
        rel_errors.append(err)

    plt.scatter(range(1, max_i+1), rel_errors)
    plt.xlabel("i - multiplicity of kmers")
    plt.ylabel("H[i] rel. errors")
    if not relative:
        plt.yscale('log')
    plt.show()


if __name__ == '__main__':
    default_params = {'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100}

    covest_estimate_on_kmerlight_and_jellyfish('error_rate', [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1],
                                               'experiment1/dataset_e{error_rate:.2f}')
    # 2% rel. error, higher relative error with for H[i] at higher i (5%)
    # 5 / 5.05 for error 5%, 4.95 / 5.20 for error 10% let us look at the histograms
    # plot_histograms('experiment1/dataset_e0.05')
    # plot_histogram_approximation_errors('experiment1/dataset_e0.05')
    #  plot_histograms('experiment1/dataset_e0.10')
    # plot_histogram_approximation_errors('experiment1/dataset_e0.10')
    # even though the histogram estimate is very accurate, large relative errors at higer i seem to influence covest

    covest_estimate_on_kmerlight_and_jellyfish('genome_size', [1000, 10000, 100000, 500000, 1000000],
                                               'experiment1/dataset_L{genome_size}')
    # 10% rel. error for higher i for small genomes
    # 5 / 5.15 for genome 1000bp long
    # plot_histograms('experiment1/dataset_L1000')
    # plot_histogram_approximation_errors('experiment1/dataset_L1000')
    # plot_histograms('experiment1/dataset_L10000')
    # plot_histogram_approximation_errors('experiment1/dataset_L10000')
    # again, the same pattern

    covest_estimate_on_kmerlight_and_jellyfish('coverage', [0.1, 0.5, 1, 2, 3, 4, 5, 10, 20],
                                               'experiment1/dataset_c{coverage}')
    # plot_histograms('experiment1/dataset_c0.1')
    # plot_histogram_approximation_errors('experiment1/dataset_c0.1') # only 3 columns
    # with long genome and relatively small error, coverage does not influence covest
    # at coverage 0.1 approximated histogram creates 1% error in estimate

    # is coverage irrelevant or is it just bad picture?
    #   0.5% - 1% estimate inaccuracy / much smaller error than when changing other parameters

    default_params['error_rate'] = 0.03
    covest_estimate_on_kmerlight_and_jellyfish('coverage', [0.1, 0.5, 1, 2, 3, 4, 5, 10, 20],
                                               'experiment1/datasetx_c{coverage}_e{error_rate:.2f}')
    default_params['error_rate'] = 0.01
    # plot_histograms('experiment1/datasetx_c0.1_e0.03') # 3 columns
    # plot_histograms('experiment1/datasetx_c0.5_e0.03') # 4 columns
    # with 3% error rate (smaller genomes are irrelevant), coverage is slightly more important
    # at coverage 0.1 20% error in estimate, at coverages0.5 - 2 about 2% error, with higher coverage less than 1%

    covest_estimate_on_kmerlight_and_jellyfish('k', [13, 17, 21, 25, 29, 33, 37, 41], 'experiment1/dataset_k{k}')
    # Covest works strange with different k-s (maybe just random error)
    # too small k=13, collisions affect covest estimate
    # covest estimate seems affected only little with different k values 

    # Further analysis

    # is really approximation error higher at the end of histogram?  (isn't it just bad picture?)
    #   yes - plot errors per histogram:
    #       relative errors increases with i
    #       (absolute errors decreases, as H[i] decreases)
    #   use relative errors plot

