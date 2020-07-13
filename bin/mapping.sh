#!/bin/bash

# ----------------------------------------------------------------------------------------------- #
# mapping.sh
# ----------------------------------------------------------------------------------------------- #
# This script will take in a sample name, a set of reference sequences, and an RNA-seq library,
# and output a mapped BAM file with the input reads mapped to the reference sequences. The script
# requires BWA and samtools, and expects 8 CPU cores. User can override this number of threads by
# changing that parameter (minimum: 1).
# ----------------------------------------------------------------------------------------------- #


# Enable strict error protection; if anything usual happens, stop & exit with an error
set -eo pipefail

# Take in the sample name, reference sequences, & input reads
sample=$1
reference=$2
reads=$3

# Set number of CPUs to use
number_of_threads=$4

# Make sure required input files are provided [if not, exit]
if [[ "$#" < 3 ]]; then
    echo "Missing required input parameters."
    echo "Usage: $0 sample-name reference.fasta reads.fastq"
    echo "Exiting with error code 1..."
    exit 1
fi

# Make sure bwa in installed [if not, exit]
command -v bwa > /dev/null || {
        echo -e "ERROR: This script requires 'bwa' but it could not found. \n" \
                "Please install this application. \n" \
                "Exiting with error code 2..." >&2; exit 2
        }

# TELL USER IS IT BEGINNING
echo "Beginning mapping of $sample to $(basename $reference '.fasta') at: $(date)"

# INDEX THE GENOME
bwa index \
-p bwa_index \
$reference

# MAP THRE READS AND REMOVE UNMAPPED READS
bwa mem \
    -t $number_of_threads \
    bwa_index \
    ${reads} | \
samtools view --threads $number_of_threads -F 4 -bh - | \
samtools sort --threads $number_of_threads - > \
    "${sample}.$(basename $reference '.fasta').mapped-and-sorted.bam"

# SAVE NUMBER OF MAPPED READS TO EACH REFERENCE SEQUENCE
samtools idxstats "${sample}.$(basename $reference '.fasta').mapped-and-sorted.bam" > \
    "${sample}.$(basename $reference '.fasta').mapping-summary.txt"

# TELL USER IT IS COMPLETE
echo "Succesfully mapped reads! Finished at $(date)"
