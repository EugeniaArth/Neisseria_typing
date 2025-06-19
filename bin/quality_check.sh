#!/bin/bash

# Define the base directory containing the country folders
BASE_DIR="/Fastq"

# Loop through each country folder
for country_dir in "$BASE_DIR"/*/; do
    # Create a directory for FastQC reports if it doesn't exist
    FASTQC_DIR="${country_dir}/fastqc_reports"
    mkdir -p "$FASTQC_DIR"

    # Find all paired fastq.gz files in the country directory
    find "${country_dir}/" -name "*_1.fastq.gz" | while read fastq_file; do
        # Run FastQC on each fastq.gz file and save the report in the fastqc_reports folder
        echo "Running FastQC on $fastq_file"
        fastqc -t 20 "$fastq_file" --outdir="$FASTQC_DIR" --threads 28
    done

    # Find all paired fastq.gz files in the country directory
    find "${country_dir}/" -name "*_2.fastq.gz" | while read fastq_file; do
        # Run FastQC on each fastq.gz file and save the report in the fastqc_reports folder
        echo "Running FastQC on $fastq_file"
        fastqc "$fastq_file" --outdir="$FASTQC_DIR" --threads 28
    done

    # Run MultiQC to aggregate FastQC reports for the current country folder
    echo "Running MultiQC in $country_dir"
    multiqc "$FASTQC_DIR" -o "$FASTQC_DIR"
done

echo "Quality control completed for all country folders."
