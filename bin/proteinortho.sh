#!/bin/bash

#dir with faa and gff3
BASE_DIR="/bakta/output"
#output eith core orthologues  .faa

OUTPUT_DIR="/proteinortho/core_files"

# Step 1 Collect only .faa files, excluding *.hypotheticals.faa
FASTA_FILES=$(find "$BASE_DIR" -mindepth 2 -type f -name "*.faa" ! -name "*.hypotheticals.faa")


# Step 2 Run proteinortho once on all of them
proteinortho -cov=90 -e=1e-10 -threads_per_process=2  -project=global_orthologs $FASTA_FILES

# Step 3 Extract core orthologues for all samples
proteinortho_grab_proteins.pl -core global_orthologs.proteinortho.tsv $FASTA_FILES

# Step 4: Move all *.faa.core files to core_files folder
find . -maxdepth 1 -type f -name "*.faa.core" -exec mv {} "$OUTPUT_DIR" \;

echo "ProteinOrtho completed. Core files moved to: $OUTPUT_DIR"
