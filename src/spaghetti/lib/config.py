import os

try:
    fastq_target_dir = os.environ['FASTQ_DIR']
    num_cpu_cores = int(os.environ['NUM_CPU_CORES'])
    reads_url = os.environ['READS_URL']
    reads_file = os.environ['READS_FILE']
    result_file = os.environ['RESULT_FILE']
except KeyError as e:
    print('Missing environment variables. Makes sure you have all set.')
    raise SystemExit
