#!/usr/bin/env bash

set -euo pipefail

PLAN="data/metadata/gse50760_download_plan_3patients.tsv"
TEMP_DIR="data/raw/test_download"
FINAL_DIR="data/raw/fastq"

mkdir -p "$TEMP_DIR"
mkdir -p "$FINAL_DIR"

tail -n +2 "$PLAN" | while IFS=$'\t' read -r sample_name srr_accession fastq_1_name fastq_2_name
do
    final_fastq_1="$FINAL_DIR/$fastq_1_name"
    final_fastq_2="$FINAL_DIR/$fastq_2_name"

    echo "----------------------------------------"
    echo "Processing sample: $sample_name"
    echo "Run accession: $srr_accession"

    if [[ -f "$final_fastq_1" && -f "$final_fastq_2" ]]; then
        echo "Final FASTQ files already exist. Skipping $sample_name."
        continue
    fi

    cd "$TEMP_DIR"

    echo "Downloading $srr_accession with prefetch..."
    prefetch "$srr_accession"

    echo "Converting $srr_accession to FASTQ..."
    fasterq-dump "$srr_accession" --split-files --threads 4

    echo "Compressing FASTQ files..."
    gzip "${srr_accession}_1.fastq"
    gzip "${srr_accession}_2.fastq"

    echo "Testing gzip integrity..."
    gzip -t "${srr_accession}_1.fastq.gz"
    gzip -t "${srr_accession}_2.fastq.gz"

    echo "Moving files to final folder..."
    mv "${srr_accession}_1.fastq.gz" "../../raw/fastq/$fastq_1_name"
    mv "${srr_accession}_2.fastq.gz" "../../raw/fastq/$fastq_2_name"

    echo "Removing temporary SRA directory..."
    rm -rf "$srr_accession"

    cd - > /dev/null

    echo "Completed: $sample_name"
	
done

echo "All planned samples processed."