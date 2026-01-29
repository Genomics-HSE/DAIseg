#!/bin/bash

# ============================================================================
# extract.samples.sh - Extract specific samples from 1000 Genomes VCF file
# ============================================================================

set -e  # Exit immediately on any error

JSON=$1
nthr=$2

# Validate input arguments
if [[ -z "$JSON" || -z "$nthr" ]]; then
    echo " Usage: $0 <config.json> <num_threads>"
    exit 1
fi

# Extract parameters from JSON
PREF="$(jq -r '.prefix' "${JSON}")"
start_time=$(date +%s)

# Get path to original 1000GP VCF
FILE_1kG_RAW="$(jq -r '.files["1000GP_files"].vcf_initial' "${JSON}")"
# Expand tilde to home directory if present
FILE_1kG="${FILE_1kG_RAW/#\~/$HOME}"

# --- CHECK 1: Verify input file exists ---
if [[ ! -f "$FILE_1kG" ]]; then
    echo "   CRITICAL ERROR: Input VCF file does not exist!"
    echo "   Path checked: $FILE_1kG"
    echo "   Please verify 'vcf_initial' in JSON or external drive mount."
    exit 1
fi

# Normalize chromosome name
CHROM_RAW=$(jq -r '.CHROM' "$JSON")
CHROM=${CHROM_RAW#chr}

# Output filenames
BASE_VCF_NAME="$(jq -r '.files["1000GP_files"].vcf' "${JSON}")"
OUTPUT_NAME=$(echo "$BASE_VCF_NAME" | sed "s/CHROMOSOME/chr${CHROM}/")
FILTERED_1kG="$PREF/$OUTPUT_NAME"

echo " Output directory: $PREF"
echo " Output VCF: $FILTERED_1kG"

# Create output directory
mkdir -p "$PREF"

echo ""
echo " Step 0: Filtering 1000 Genomes samples..."

# Create temporary sample list
SAMPLE_LIST=$(mktemp)
jq -r '.samples.ingroup[], .samples.outgroup[]' "$JSON" > "$SAMPLE_LIST"

SAMPLE_COUNT=$(wc -l < "$SAMPLE_LIST")
echo " Extracting $SAMPLE_COUNT samples..."

if [[ "$SAMPLE_COUNT" -eq 0 ]]; then
    echo " Error: Sample list is empty. Check JSON configuration."
    rm "$SAMPLE_LIST"
    exit 1
fi

# --- RUN BCFTOOLS to extract samples ---
echo "  Running bcftools..."
bcftools view --threads "$nthr" -S "$SAMPLE_LIST"  "$FILE_1kG" -Oz -o "$FILTERED_1kG"

# --- CHECK 2: Verify output was created successfully ---
if [[ ! -s "$FILTERED_1kG" ]]; then
    echo " CRITICAL ERROR: Output VCF is empty or failed to write."
    exit 1
fi

echo " Indexing VCF..."
bcftools index --tbi --threads "$nthr" "$FILTERED_1kG"

# Cleanup temporary files
rm "$SAMPLE_LIST"

end_time=$(date +%s)
duration=$((end_time - start_time))
echo ""
echo " Success! Filtered VCF saved to: $FILTERED_1kG"
echo "  Total time: ${duration} seconds"
