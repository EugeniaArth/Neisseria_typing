#!/bin/bash

date  ## echo the date at start

export PATH=$HOME/sratoolkit.3.0.7-mac64/bin:$PATH

#  output directory for Fastq files
BASE_OUTPUT_DIR="/Fastq"
LOG_FILE="download_log.txt"
SUMMARY_FILE="/Fastq/summary.txt"

# Create base output directory if it doesn't exist
mkdir -p "$BASE_OUTPUT_DIR"

# Initialize summary file with headers
echo -e "Country\tRun Accession\tResult" > $SUMMARY_FILE
# Create or clear the log file
> $LOG_FILE

# Skip the first line of the input file (header) and read subsequent lines
tail -n +2 for_download.txt | while IFS=$'\t' read -r country run_accession; do
    # Remove any trailing carriage return or other invisible characters
    country=$(echo "$country" | tr -d '\r')
    run_accession=$(echo "$run_accession" | tr -d '\r')

    # Define directories for the country and run accession
    COUNTRY_DIR="$BASE_OUTPUT_DIR/$country"

    # Create directories for the country if they don't exist
    mkdir -p "$COUNTRY_DIR"

    # Log start of download for current accession
    echo "Processing $run_accession for $country..." | tee -a $LOG_FILE

    # Remove any existing lock file before starting a new download
    LOCK_FILE="$COUNTRY_DIR/$run_accession.sra.lock"
    if [ -f "$LOCK_FILE" ]; then
        echo "Removing existing lock file: $LOCK_FILE" | tee -a $LOG_FILE
        rm -f "$LOCK_FILE"
    fi

    # Attempt to prefetch the SRA file directly into the COUNTRY_DIR without nested folder
    prefetch $run_accession -O "$COUNTRY_DIR" 2>> $LOG_FILE
    if [[ $? -ne 0 ]]; then
        # Log failure in summary and skip this run
        echo "$country,$run_accession,Failed (prefetch)" >> $SUMMARY_FILE
        continue
    fi

    # Change to the country directory (not the nested one)
    cd "$COUNTRY_DIR"

    # Attempt to convert the SRA file to FASTQ
    fasterq-dump "$run_accession" --threads 20 --split-files 2>> "$LOG_FILE"
    if [[ $? -ne 0 ]]; then
        # Log failure in summary and skip this run
        echo "$country,$run_accession,Failed (fasterq-dump)" >> "$SUMMARY_FILE"
        cd "$BASE_OUTPUT_DIR"
        continue
    fi

    # Compress the FASTQ files
    gzip -f ${COUNTRY_DIR}/${run_accession}_1.fastq 2>> "$LOG_FILE"
    gzip -f ${COUNTRY_DIR}/${run_accession}_2.fastq 2>> "$LOG_FILE"


    # Remove temporary and unnecessary files
    rm -r "$COUNTRY_DIR/$run_accession" 2>> "$LOG_FILE"

    # Log success in the summary
    echo "$country,$run_accession,Success" >> "$SUMMARY_FILE"

    # Return to the base directory
    cd "$BASE_OUTPUT_DIR"

done

echo "Download process completed. Check $SUMMARY_FILE and $LOG_FILE for details."
