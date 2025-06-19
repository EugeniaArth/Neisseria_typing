#!/bin/bash

# Export the path to BWA
export PATH="/Users/eugenianikonorova/bwa:$PATH"
#create indexes for reference genomes
#bwa index -p FA19 /Users/eugenianikonorova/Genomes_mol_typing/Neisseria/reference/FA19.fasta


# Path to the reference genome index
REFERENCE_INDEX="/reference_folder/FA19"

# Base directories
FASTQ_DIR="/Fastq"
MAPPED_BASE_DIR="/GATK4/mapped/"
LOG_DIR="/GATK4/logs"
RESULTS_FILE="/GATK4/mapping_results.txt"

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"
mkdir -p "$MAPPED_BASE_DIR"

# Initialize results file with header
echo -e "Sample Name\tCountry\tMapped" > "$RESULTS_FILE"

# Iterate over each country directory
for COUNTRY_DIR in $FASTQ_DIR/*; do
    # Ensure it is a directory
    if [ -d "$COUNTRY_DIR" ]; then
        COUNTRY_NAME=$(basename "$COUNTRY_DIR")

        # Create the mapped directory for the country if it doesn't exist
        MAPPED_COUNTRY_DIR="$MAPPED_BASE_DIR/$COUNTRY_NAME"
        mkdir -p "$MAPPED_COUNTRY_DIR"

        # Iterate over each fastq.gz file pair in the Trimmed subfolder
        for F in "$COUNTRY_DIR/"*_1.fastq.gz; do
            I=${F%%_1.fastq.gz}
            SAMPLE_NAME=$(basename "$I")
            LOG_FILE="$LOG_DIR/${SAMPLE_NAME}_mapping.log"

            # Perform bwa mem alignment and log the output
            #bwa mem -t  30 "$REFERENCE_INDEX" "${I}_trimmed_1.fastq.gz" "${I}_trimmed_2.fastq.gz" > "${I}.sam" 2> "$LOG_FILE"
            bwa mem -t  30 "$REFERENCE_INDEX" "${I}_1.fastq.gz" "${I}_2.fastq.gz" > "${I}.sam" 2> "$LOG_FILE"
            # Check if mapping was successful by inspecting the exit status
            if [ $? -eq 0 ]; then
                MAPPED_STATUS="yes"
                echo "${SAMPLE_NAME} mapping is finished" >> "$LOG_FILE"

                # Convert SAM to BAM
                samtools view -bS -F 0x4 "${I}.sam" -o "${I}.bam" >> "$LOG_FILE" 2>&1
                echo "${SAMPLE_NAME}.bam is received" >> "$LOG_FILE"

                # Sort BAM file
                samtools sort "${I}.bam" -o "${I}.sorted.bam" >> "$LOG_FILE" 2>&1
                echo "${SAMPLE_NAME}.bam is sorted" >> "$LOG_FILE"

                # Move the BAM files to the mapped directory
                mv "${I}.bam" "${MAPPED_COUNTRY_DIR}/"
                mv "${I}.sorted.bam" "${MAPPED_COUNTRY_DIR}/"
                echo "${SAMPLE_NAME} was moved to $MAPPED_COUNTRY_DIR" >> "$LOG_FILE"
            else
                MAPPED_STATUS="no"
                echo "${SAMPLE_NAME} mapping failed" >> "$LOG_FILE"
            fi

            # Clean up SAM file
            rm "${I}.sam"
            echo "${SAMPLE_NAME}.sam was removed" >> "$LOG_FILE"

            # Save the result to the results file
            echo -e "${SAMPLE_NAME}\t${COUNTRY_NAME}\t${MAPPED_STATUS}" >> "$RESULTS_FILE"
        done
    fi
done
