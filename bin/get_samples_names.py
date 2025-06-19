import os

#grab the sample_name (from sample.fasta) from folders to next step

# Define the directory containing the country folders
root_directory = 'dir_with_countries_folders_and_samples'  #base folder with fasta and fas
output_file = 'fasta_downloaded.txt'
folders_to_ignore = ['All']  # List of folders to ignore

with open(output_file, 'w') as outfile:
    for country_folder in os.listdir(root_directory):
        # Skip folders that should be ignored
        if country_folder in folders_to_ignore:
            continue

        country_path = os.path.join(root_directory, country_folder)
        if os.path.isdir(country_path):
            for file_name in os.listdir(country_path):
                if (file_name.endswith('.fas') or file_name.endswith('.fasta')):
                    sample_name = os.path.splitext(file_name)[0]
                    if '_' in sample_name:
                        sample_name = sample_name.split('_')[0]
                    outfile.write(f"{country_folder}\t{sample_name}\n")

print(f"Sample names and country names have been written to {output_file}")
