from collections import defaultdict
import os
import resource
import subprocess

import numpy as np

covest_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'covest')
kmerlight_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'kmerlight')

sequence_generator = os.path.join(covest_path, 'tools', 'simulator', 'generate_sequence.py') +\
                     ' -g {genome_size} {base_dir}/genome.fa'
reads_generator = os.path.join(covest_path, 'tools', 'simulator', 'read_simulator.py') +\
                  ' -c {coverage} -e {error_rate} -r {read_length} {base_dir}/genome.fa {base_dir}/reads.fa'

jellyfish_count = 'jellyfish count -m {k} -s 100M -t 4 -C {base_dir}/reads.fa -o {base_dir}/table.jf'
jellyfish_hist = 'jellyfish histo {base_dir}/table.jf -o {base_dir}/jellyfish.dist'
run_kmerlight = os.path.join(kmerlight_path, 'kmerlight') +\
                ' -k {k} -f {kmerlight_f} -a {kmerlight_alternative}' \
                ' -o {base_dir}/kmerlight_partial.dist -i {reads_dir}/reads.fa'

run_covest = 'covest {base_dir}/{histogram} -rs {read_length} -k {k} -r {read_length} -m basic -sf 1'
run_covest_alternative = 'covest {base_dir}/{histogram} -rs {read_length} -k {k} -r {read_length} -m basic' \
                         ' -lm {likelihood} -lf {likelihood-data-file} -sf 1'


def expected_number_of_kmers(coverage, genome_size, read_size, k):
    return int(genome_size / read_size * coverage) * (read_size - k + 1)


def run(command, output=None, errors=None):
    f = open(output, 'w') if output is not None else None
    g = open(errors, 'w') if errors is not None else None
    return subprocess.call(command.split(), stdout=f, stderr=g)


_utime = 0
_stime = 0


def start_stopwatch():
    info = resource.getrusage(resource.RUSAGE_CHILDREN)
    global _utime, _stime
    _utime = info.ru_utime
    _stime = info.ru_stime


def stop_stopwatch(label):
    info = resource.getrusage(resource.RUSAGE_CHILDREN)
    user = info.ru_utime - _utime
    system = info.ru_stime - _stime
    print('{} user: {:.3f}, system: {:.3f}'.format(label, user, system))


def generate_reads(params):
    run(sequence_generator.format(**params))
    run(reads_generator.format(**params))


def normalize_kmerlight_output(infile, outfile):
    hist = []
    with open(infile, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            try:
                l = line.split()
                i = int(l[0].strip('f'))
                cnt = int(l[1])
                if cnt == 0:
                    continue
                hist.append((i, cnt))
            except:
                pass

    with open(outfile, 'w') as f:
        for i, cnt in hist:
            f.write('{} {}\n'.format(i, cnt))


def calculate_jellyfish_histogram(params):
    run(jellyfish_count.format(**params))
    run(jellyfish_hist.format(**params))


def calculate_kmerlight_histogram(params, outfile='kmerlight'):
    if 'kmerlight_f' not in params:
        params['kmerlight_f'] = 100

    if 'reads_dir' not in params:
        params['reads_dir'] = params['base_dir']

    if 'kmerlight_alternative' not in params:
        params['kmerlight_alternative'] = 0

    run(run_kmerlight.format(**params),
        output='{base_dir}/{outfile}.log'.format(**params, outfile=outfile))
    hist_path = '{base_dir}/kmerlight_partial.dist'.format(**params)
    out_path = '{base_dir}/{outfile}.dist'.format(**params, outfile=outfile)
    normalize_kmerlight_output(hist_path, out_path)


def load_histogram(infile):
    hist = defaultdict(int)
    with open(infile, 'r') as f:
        for line in f:
            l = line.split()
            i = int(l[0])
            cnt = int(l[1])
            if cnt > 0:
                hist[i] = cnt
    return hist


def load_histograms_to_matrix(filenames, maxf=None):
    if maxf is None:
        maxf = max(load_histogram(filenames[0]).keys())

    histograms = np.zeros([len(filenames), maxf + 1])

    for i, filename in enumerate(filenames):
        hist = load_histogram(filename)
        for k, v in hist.items():
            if k <= maxf:
                histograms[i, k] = v

    return histograms


def run_covest_on_histograms(params, identifiers):
    mparams = dict(**params)

    for name in identifiers:
        mparams['histogram'] = name + '.dist'
        run(run_covest.format(**mparams),
            output='{base_dir}/{name}.est'.format(**mparams, name=name),
            errors='{base_dir}/{name}.est.log'.format(**mparams, name=name))


def run_covest_on_all_histograms(params):
    mparams = dict(**params)

    for file in os.listdir(params['base_dir']):
        if file.endswith('.dist'):
            mparams['histogram'] = file
            name, ext = os.path.splitext(file)
            run(run_covest.format(**mparams),
                output='{base_dir}/{filename}.est'.format(**params, filename=name),
                errors='{base_dir}/{filename}.est.log'.format(**params, filename=name))


def load_covest_results(infile):
    results = dict()
    with open(infile, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            k, v = line.split(':')
            results[k.strip()] = v.strip()
    return results


def load_coverage_estimates(filenames):
    estimates = []
    for filename in filenames:
        estimates.append(float(load_covest_results(filename).get('coverage', 0)))
    return np.array(estimates)


def save_error_statistics(dirname, mean_errors, sds):
    with open('{}/error_stats.errs'.format(dirname), 'w') as f:
        f.write(' '.join(map(str, mean_errors)))
        f.write('\n')
        f.write(' '.join(map(str, sds)))
        f.write('\n')


def load_error_statistics(dirname):
    with open('{}/error_stats.errs'.format(dirname), 'r') as f:
        mean_errors = [float(x) for x in f.readline().split()]
        sds = [float(x) for x in f.readline().split()]
        return np.array(mean_errors), np.array(sds)
