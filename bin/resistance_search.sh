#!/bin/bash

INPUT_DIR="Fasta_analyzed"
OUTPUT_DIR="/rgi_results"
LOG_FILE="${OUTPUT_DIR}/rgi_run.log"

mkdir -p "$OUTPUT_DIR"
echo "RGI run started: $(date)" > "$LOG_FILE"

find "$INPUT_DIR" -type f \( -iname "*.fasta" -o -iname "*.fas" -o -iname "*.faa" \) | while read -r file; do
    country=$(basename "$(dirname "$file")")
    sample=$(basename "$file" | sed 's/\.[^.]*$//')
    out_dir="${OUTPUT_DIR}/${country}"
    mkdir -p "$out_dir"

    echo "Processing $sample from $country"
    echo "[$(date)] Processing $sample from $country" >> "$LOG_FILE"

    docker run --rm \
      -v "$(dirname "$file")":/in \
      -v "$out_dir":/out \
      quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0 \
      rgi main \
      --input_sequence "/in/$(basename "$file")" \
      --output_file "/out/${sample}_rgi" \
      --input_type contig \
      > "${out_dir}/${sample}_rgi.log" 2>&1

    STATUS=$?

    if [ $STATUS -eq 0 ]; then
        echo "✅ $sample completed successfully."
        echo "[$(date)] ✅ $sample completed successfully." >> "$LOG_FILE"
        echo "Output files for $sample:"

        # 🔥 Delete everything except final results and log
        find "$out_dir" -type f -name "${sample}*.temp*" -delete
        find "$out_dir" -type f -name "${sample}*.blastRes.*" -delete
        find "$out_dir" -type f -name "${sample}*.db.*" -delete
    else
        echo "❌ $sample failed (exit code $STATUS). See log: ${out_dir}/${sample}_rgi.log"
        echo "[$(date)] ❌ $sample failed (exit code $STATUS)" >> "$LOG_FILE"
    fi
done

echo "RGI run finished: $(date)" >> "$LOG_FILE"
