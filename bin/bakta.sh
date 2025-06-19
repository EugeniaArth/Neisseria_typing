#!/bin/bash
shopt -s nullglob extglob

export PATH="/opt/homebrew/Cellar/aragorn/1.2.41/bin:$PATH"

# Input parent folder (contains country subfolders)
FASTA_PARENT="/Fasta_analyzed"

# Bakta DB path
DB="/bakta/db-light"

# Output folder
OUTPUT="/bakta/output"

# Loop through all country folders
for COUNTRY_DIR in "$FASTA_PARENT"/*; do
    if [[ -d "$COUNTRY_DIR" ]]; then
        COUNTRY=$(basename "$COUNTRY_DIR")
        echo "📂 Processing country: $COUNTRY"

        # Find all .fas and .fasta files in this country folder
        for fasta_file in "$COUNTRY_DIR"/*.@(fasta|fas); do
            [[ -e "$fasta_file" ]] || continue  # Skip if no files

            base_name=$(basename "$fasta_file")
            base_name="${base_name%.*}"

            # Define output path for this sample
            SAMPLE_OUT="${OUTPUT}/${COUNTRY}/${base_name}"
            mkdir -p "$SAMPLE_OUT"

            echo "🔬 Annotating $base_name..."

            bakta --db "$DB" --force \
                  --output "$SAMPLE_OUT" \
                  --threads 8 \
                  --skip-crispr \
                  --skip-plot \
                  --prefix "$base_name" \
                  "$fasta_file"
        done
    fi
done

echo "🎉 All annotations and extractions complete!"
