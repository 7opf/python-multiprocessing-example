from multiprocessing import Pool
from lib.config import num_cpu_cores, reads_url, reads_file, result_file
from lib.experiment import Experiment
import urllib
import pandas as pd
import os.path
from filelock import FileLock


def do_experiment(arg_tuple):
    # run the experiment
    exp = Experiment(arg_tuple.run_accession, arg_tuple.fastq_ftp, arg_tuple.fastq_md5)
    # acquire a lock on the result file (needed since this function is concurrently executed)
    lock = FileLock(result_file + '.lock')
    with lock:
        # append the experiment result as a row to the tab-separated values file
        row = exp.run_accession + '\t' + exp.fastq_dl_success + '\t' + exp.cortex_conversion_success + '\t' + \
              exp.kraken_file_success + '\t' + exp.brack_file_success + '\n'
        open(result_file, "a").write(row)


def load_reads_file():
    # download the whole reads file to disk from the given URL
    # to the given path if it isn't already there
    if not os.path.isfile(reads_file):
        print 'DOWNLOADING READS FILE...'
        urllib.urlretrieve(reads_url, reads_file)
    # load the CSV into a DataFrame (tab delimited, utf-8 line breaks)
    print 'LOADING READS FILE FROM DISK...'
    return pd.read_csv(reads_file, sep='\t', lineterminator='\n')


if __name__ == '__main__':
    print 'STARTING MULTI-PROCESS JOB ACROSS', num_cpu_cores, 'CORES'

    # init the reads_table dataframe
    reads_table = load_reads_file()

    # init the result file (overwrites!)
    headers = 'run_accession\tfastq_dl_success\tcortex_conversion_success\tkraken_file_success\tbrack_file_success\n'
    with open(result_file, 'w') as f:
        f.write(headers)

    # init the process pool to which experiments are sent to
    # https://docs.python.org/2/library/multiprocessing.html
    # https://www.ellicium.com/python-multiprocessing-pool-process/
    process_pool = Pool(num_cpu_cores)

    # map the data frame iterable to the experiment function.
    # in other words; create a task out of each row of the table
    # by passing the row as arguments to the experiment function.
    # the multiprocessing library deals with scheduling the tasks
    # between the cores
    process_pool.map(do_experiment, reads_table.itertuples())

    print 'COMPLETED MULTI-PROCESS JOB'
