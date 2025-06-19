import os
import pandas as pd

# Base directory containing country folders
BASE_DIR = "/Fastq"

# Output file for failed samples
output_file = os.path.join(BASE_DIR, "failed_samples.txt")

# Columns to check for "fail"
columns_to_check = ["per_base_sequence_quality", "per_sequence_gc_content", "adapter_content"]

# Prepare output list for failed samples
failed_samples = []

# Iterate through each country folder
for country_folder in os.listdir(BASE_DIR):
    country_dir = os.path.join(BASE_DIR, country_folder)
    fastqc_reports_dir = os.path.join(country_dir, "fastqc_reports")
    multiqc_file = os.path.join(fastqc_reports_dir, "multiqc_data", "multiqc_fastqc.txt")

    # Check if multiqc file exists
    if os.path.exists(multiqc_file):
        # Read the multiqc_fastqc.txt file
        df = pd.read_csv(multiqc_file, sep="\t")

        # Iterate through the rows of the dataframe
        for index, row in df.iterrows():
            # Check if any of the specified columns have "fail"
            if any(row[col] == 'fail' for col in columns_to_check):
                # Append the sample, country, and results of each column to the list
                failed_samples.append([
                    row['Sample'],
                    country_folder,
                    row['per_base_sequence_quality'],
                    row['per_sequence_gc_content'],
                    row['adapter_content']
                ])

# Create DataFrame from failed samples
failed_df = pd.DataFrame(failed_samples, columns=['sample', 'country', 'per_base_sequence_quality', 'per_sequence_gc_content', 'adapter_content'])

# Save to a text file
failed_df.to_csv(output_file, sep="\t", index=False)

print(f"Failed samples report saved to {output_file}")
