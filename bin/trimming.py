import os
import subprocess

# Base directory containing country folders
BASE_DIR = "/Fastq"

# Directory to save trimmed files
TRIMMED_DIR = "/Trimmed"
RESULTS_FILE = "trimming_results.txt"

# Initialize results file with header
if not os.path.exists(RESULTS_FILE):
    with open(RESULTS_FILE, "w") as results:
        results.write("Sample Name\tCountry\tStatus\n")

# Create trimming directory if it doesn't exist
if not os.path.exists(TRIMMED_DIR):
    os.makedirs(TRIMMED_DIR)

# Adapter file path for Trimmomatic
ADAPTERS_PATH = "Sequencing_adaptors.fasta"

# Trimmomatic command template
trimmomatic_cmd_template = (
    "trimmomatic PE -threads 20 -phred33 {input_1} {input_2} "
    "{output_1_paired} {output_1_unpaired} {output_2_paired} {output_2_unpaired} "
    "ILLUMINACLIP:{adapters}:2:30:10 LEADING:25 TRAILING:25 SLIDINGWINDOW:10:30 MINLEN:30"
)

# Iterate through country directories
for country in os.listdir(BASE_DIR):
    country_dir = os.path.join(BASE_DIR, country)
    if not os.path.isdir(country_dir):
        continue

    trimming_dir = os.path.join(TRIMMED_DIR, country)

    # Create trimming directory if it doesn't exist
    if not os.path.exists(trimming_dir):
        os.makedirs(trimming_dir)

    # Iterate through files in the country directory
    for file_name in os.listdir(country_dir):
        if not file_name.endswith("_1.fastq.gz"):
            continue

        sample_prefix = file_name.rsplit("_", 1)[0]
        input_1 = os.path.join(country_dir, f"{sample_prefix}_1.fastq.gz")
        input_2 = os.path.join(country_dir, f"{sample_prefix}_2.fastq.gz")

        output_1_paired = os.path.join(trimming_dir, f"{sample_prefix}_trimmed_1.fastq.gz")
        output_1_unpaired = os.path.join(trimming_dir, f"{sample_prefix}_unpaired_1.fastq.gz")
        output_2_paired = os.path.join(trimming_dir, f"{sample_prefix}_trimmed_2.fastq.gz")
        output_2_unpaired = os.path.join(trimming_dir, f"{sample_prefix}_unpaired_2.fastq.gz")

        log_file = os.path.join(trimming_dir, f"{sample_prefix}_trimming.log")

        trimmomatic_cmd = trimmomatic_cmd_template.format(
            input_1=input_1,
            input_2=input_2,
            output_1_paired=output_1_paired,
            output_1_unpaired=output_1_unpaired,
            output_2_paired=output_2_paired,
            output_2_unpaired=output_2_unpaired,
            adapters=ADAPTERS_PATH
        )

        print(f"Processing Country: {country}, Sample: {sample_prefix}")
        try:
            with open(log_file, "w") as log:
                result = subprocess.run(trimmomatic_cmd, shell=True, stdout=log, stderr=log)

            if result.returncode == 0:
                status = "Trimmed"
                print(f"Country: {country}, Sample: {sample_prefix} - Completed successfully")
            else:
                status = "Failed"
                print(f"Country: {country}, Sample: {sample_prefix} - Failed")

        except Exception as e:
            status = "Failed"
            print(f"Country: {country}, Sample: {sample_prefix} - Failed")
            print("Exception occurred:", e)

        with open(RESULTS_FILE, "a") as results:
            results.write(f"{sample_prefix}\t{country}\t{status}\n")

print("Trimming completed for all samples.")
