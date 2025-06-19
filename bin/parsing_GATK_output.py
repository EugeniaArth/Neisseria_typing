import pandas as pd
import re
import os

# Function to extract field values from the info column
def extract_field(info, field_name):
    pattern = rf'{field_name}\s*"([^"]+)"'
    match = re.search(pattern, info)
    return match.group(1) if match else None

def extract_gene_name(info):
    if "gbkey=Gene" in info:
        match = re.search(r"Name=([^;]+)", info)
        return match.group(1) if match else None

def extract_gene_biotype(info):
    if "gbkey=Gene" in info:
        match = re.search(r"gene_biotype=([^;]+)", info)
        return match.group(1) if match else None

def extract_protein_name(info):
    if "gbkey=CDS" in info:
        match = re.search(r"product=([^;]+)", info)
        return match.group(1) if match else None

def extract_protein_id(info):
    if "gbkey=CDS" in info:
        match = re.search(r"ID=cds-([^;]+)", info)
        return match.group(1) if match else None

def extract_RNA(info):
    if "gbkey=rRNA" in info:
        match = re.search(r"locus_tag=([^;]+)", info)
        return match.group(1) if match else None

def extract_locus_tag(info):
    match = re.search(r"locus_tag=([^;]+)", info)
    return match.group(1) if match else None

def extract_SNPg(SNPg):
    pattern = r"c\.([A-Z]\d+[A-Z])"
    match = re.search(pattern, SNPg)
    return match.group(1) if match else None

def extract_SNPp(SNPp):
    pattern = r"p\.([A-Z]\d+[A-Z])"
    match = re.search(pattern, SNPp)
    return match.group(1) if match else None


# Process a single SNP file
def process_exonic_variant_function(file_path, output_dir):
    print(f"Processing file: {file_path}")
    filename = os.path.basename(file_path)
    prefix = '.'.join(filename.split('.')[:2])

    if os.path.getsize(file_path) == 0:
        print(f"Skipping empty file: {file_path}")
        return

    try:
        data = pd.read_csv(file_path, sep='\t', header=None)
        print(f"File read successfully: {file_path}, rows: {len(data)}")
    except pd.errors.EmptyDataError:
        print(f"Empty or malformed data in file: {file_path}")
        return
    except Exception as e:
        print(f"Error reading file: {file_path}, Error: {e}")
        return

    gene_info_dict = {}
    processed_data = []

    for _, row in data.iterrows():
        info = str(row.iloc[29])  # Column 29 contains detailed info
        snp_type = row.iloc[1]
        position = row.iloc[5]
        SNPg = row.iloc[2]
        SNPp = row.iloc[2]
        start = row.iloc[24]
        end = row.iloc[25]

        if 'gbkey=Gene' in info:
            gene_name = extract_gene_name(info) or "Unknown gene"
            gene_id = extract_locus_tag(info) or "Unknown_gene_id"
            gene_biotype = extract_gene_biotype(info) or "None"

            if gene_id:
                gene_info_dict[gene_id] = {
                    'gene_name': gene_name,
                    'gene_biotype': gene_biotype,
                    'start': start,
                    'end': end
                }
        elif 'gbkey=CDS' in info:
            rRNA = extract_RNA(info) or "None"

        elif 'gbkey=CDS' in info:
            protein_name = extract_protein_name(info) or "None"
            protein_id = extract_protein_id(info) or "None"
            locus_tag = extract_locus_tag(info)
            SNPg_name = extract_SNPg(SNPg) or "None"
            SNPp_name = extract_SNPp(SNPp) or "None"

            gene_id = extract_locus_tag(info)
            if gene_id in gene_info_dict:
                gene_entry = gene_info_dict[gene_id]
                processed_data.append([
                    gene_entry['start'], gene_entry['end'], locus_tag, gene_entry['gene_name'],
                    gene_entry['gene_biotype'], protein_name, protein_id, snp_type,
                    position, SNPg_name, SNPp_name
                ])
            else:
                processed_data.append([
                    start, end, locus_tag, gene_id, "None", protein_name, protein_id, snp_type,
                    position, SNPg_name, SNPp_name
                ])

    columns = ["Start", "End", "Locus", "Gene name", "Gene biotype", "Protein name", "Protein ID",
               "Type of SNP", "Position", "SNP", "AA change"]
    result_df = pd.DataFrame(processed_data, columns=columns)

    output_path = os.path.join(output_dir, f"{prefix}_processed_data.csv")
    os.makedirs(output_dir, exist_ok=True)
    result_df.to_csv(output_path, index=False)
    print(f"Processed data saved to: {output_path}")


# Function to iterate through multiple folders
def process_multiple_folders(base_dir):
    # Define the root output directory
    root_output_dir = '/GATK4/analysis/'

    for country_folder in os.listdir(base_dir):
        country_path = os.path.join(base_dir, country_folder)
        if os.path.isdir(country_path):
            vcf_annotated_folder = os.path.join(country_path, 'vcf_annotated')
            if os.path.exists(vcf_annotated_folder):
                output_dir = os.path.join(root_output_dir, country_folder)
                for file in os.listdir(vcf_annotated_folder):
                    if file.endswith('.snp.refGene.exonic_variant_function'):
                        file_path = os.path.join(vcf_annotated_folder, file)
                        process_exonic_variant_function(file_path, output_dir)

# Example usage
base_dir = '/GATK4/calling/'  # Adjust this path as needed
process_multiple_folders(base_dir)
