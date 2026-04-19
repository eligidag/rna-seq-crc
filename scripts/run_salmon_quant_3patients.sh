#!/usr/bin/env bash

set -euo pipefail

PLAN="data/metadata/gse50760_download_plan_3patients.tsv"
INDEX="reference/gencode_v49_salmon/salmon_index"
FASTQ_DIR="data/raw/fastq"
OUTDIR="results/salmon_quant_3patients"
THREADS=4

mkdir -p "$OUTDIR"

tail -n +2 "$PLAN" | while IFS=$'\t' read -r sample_name srr_accession fastq_1_name fastq_2_name
do
    fastq_1="$FASTQ_DIR/$fastq_1_name"
    fastq_2="$FASTQ_DIR/$fastq_2_name"
    sample_outdir="$OUTDIR/$sample_name"
    quant_file="$sample_outdir/quant.sf"

    echo "----------------------------------------"
    echo "Sample: $sample_name"
    echo "FASTQ 1: $fastq_1"
    echo "FASTQ 2: $fastq_2"

    if [[ -f "$quant_file" ]]; then
        echo "quant.sf already exists. Skipping $sample_name."
        continue
    fi

    if [[ ! -f "$fastq_1" || ! -f "$fastq_2" ]]; then
        echo "Missing FASTQ file(s) for $sample_name"
        exit 1
    fi

    echo "Running salmon quant for $sample_name ..."
    salmon quant \
        -i "$INDEX" \
        -l A \
        -1 "$fastq_1" \
        -2 "$fastq_2" \
        -p "$THREADS" \
        -o "$sample_outdir"

    echo "Completed: $sample_name"
done

echo "All planned samples checked/processed."