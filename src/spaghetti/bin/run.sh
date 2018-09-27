#!/usr/bin/env bash

# TODO: Ed - youll need to change these variables

# directory where fastq gz files will be downloaded into
FASTQ_DIR="C:\code\spaghetti\test"

# number of CPU cores to use
NUM_CPU_CORES=2

# URL of the reads text file (can be the whole txt file)
READS_URL="https://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_tree(2)%20AND%20first_created%3E=2017-01-01%22&limit=358768&length=358768&offset=1&display=report&result=read_run&fields=fastq_ftp,fastq_md5&download=txt",

# path to where the reads file will be saved (including filename)
READS_FILE="C:\code\spaghetti\test\example.ena.txt"

# path to where the result file will be saved (including filename)
RESULT_FILE="C:\code\spaghetti\test\result.txt"

# run the python main file (the variables above are loaded by src/lib/config.py)
FASTQ_DIR=${FASTQ_DIR} NUM_CPU_CORES=${NUM_CPU_CORES} READS_URL=${READS_URL} READS_FILE=${READS_FILE} RESULT_FILE=${RESULT_FILE} python ../__main__.py