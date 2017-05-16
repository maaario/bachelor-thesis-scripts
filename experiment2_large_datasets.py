# Experiment 2, generate large datasets - to be run on grid engine

from experiment2 import generate_dataset

if __name__ == '__main__':
    param_sets = [
        # {'genome_size': 1000000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        {'genome_size': 10000000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        # {'genome_size': 100000000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        #
        # {'genome_size': 100000000, 'coverage': 0.5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        # {'genome_size': 100000000, 'coverage': 2, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        # {'genome_size': 100000000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        #
        # {'genome_size': 100000000, 'coverage': 5, 'error_rate': 0.01, 'k': 21, 'read_length': 100},
        # {'genome_size': 100000000, 'coverage': 5, 'error_rate': 0.05, 'k': 21, 'read_length': 100},
        # {'genome_size': 100000000, 'coverage': 5, 'error_rate': 0.1, 'k': 21, 'read_length': 100},
    ]

    for params in param_sets:
        dirname_temp = 'experiment2/dataset_L{genome_size}_c{coverage}_e{error_rate}_k{k}'
        trial_dirname_temp = dirname_temp + '/trial{trial:02d}'
        params['kmerlight_f'] = 100
        generate_dataset(params, trial_dirname_temp, 20)
