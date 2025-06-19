import pandas as pd
import os
from collections import defaultdict

def process_summary(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for country_folder in os.listdir(input_dir):
        country_path = os.path.join(input_dir, country_folder)
        if os.path.isdir(country_path):
            summary_data = []

            for file in os.listdir(country_path):
                if file.endswith('_processed_data.csv'):
                    sample_name = file.split('.')[0]
                    file_path = os.path.join(country_path, file)
                    data = pd.read_csv(file_path)

                    grouped = data.groupby([
                        "Start", "End", "Locus", "Gene name", "Gene biotype", "Protein name", "Protein ID",
                        "Type of SNP"
                    ]).size()

                    # Create a dictionary to track counts per type of SNP
                    count_tracker = defaultdict(lambda: defaultdict(int))

                    for index, count in grouped.items():
                        start, end, locus, gene_name, gene_biotype, protein_name, protein_id, snp_type = index

                        if gene_biotype == 'unknown':
                            gene_biotype = 'Unknown'

                        count_tracker[(start, end, locus, gene_name, gene_biotype, protein_name, protein_id)][snp_type] += count

                    for key, counts in count_tracker.items():
                        start, end, locus, gene_name, gene_biotype, protein_name, protein_id = key
                        nonsynonymous_count = counts.get('nonsynonymous SNV', 0)
                        synonymous_count = counts.get('synonymous SNV', 0)

                        if nonsynonymous_count == 0 and synonymous_count == 0:
                            ratio = "not applicable, single SNP"
                        elif synonymous_count == 0:
                            ratio = "not applicable, single SNP"
                        else:
                            ratio = round(nonsynonymous_count / synonymous_count, 2)

                        summary_data.append([
                            sample_name, gene_name, locus, gene_biotype, protein_name, protein_id, start, end, ratio,
                            nonsynonymous_count, synonymous_count
                        ])

            if summary_data:
                summary_df = pd.DataFrame(summary_data, columns=[
                    "Sample Name", "Gene name", "Locus", "Gene biotype", "Protein name", "Protein ID", "Start", "End",
                    "Ratio", "Total nonsynSNP", "Total synSNP"
                ])

                output_file = os.path.join(output_dir, f"{country_folder}_summary_table.csv")
                summary_df.to_csv(output_file, index=False)
                print(f"Summary saved to {output_file}")

def combine_summaries(base_dir, output_dir):
    combined_data = []
    os.makedirs(output_dir, exist_ok=True)

    for summary_file in os.listdir(base_dir):
        if summary_file.endswith("_summary_table.csv"):
            file_path = os.path.join(base_dir, summary_file)
            print(f"Found summary file: {file_path}")
            summary_df = pd.read_csv(file_path)
            country_name = summary_file.split("_summary_table")[0]
            summary_df['Country'] = country_name
            combined_data.append(summary_df)

    if combined_data:
        combined_df = pd.concat(combined_data)
        global_output_file = os.path.join(output_dir, "combined_summary.csv")
        combined_df.to_csv(global_output_file, index=False)
        print(f"Combined summary saved to {global_output_file}")
    else:
        print("Error: No summary files were found to combine.")

#  usage
base_dir = '/GATK4/analysis/'
summary_output_dir = '/GATK4/summary/'  # Replace with the desired output folder for summaries
combine_output_dir = summary_output_dir  # Combine summaries in the same output directory

# Generate summaries from existing processed data
process_summary(base_dir, summary_output_dir)

# Combine summaries into a global file
combine_summaries(summary_output_dir, combine_output_dir)
