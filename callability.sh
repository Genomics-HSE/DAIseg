#!/bin/bash

# ============================================================================
# callability.sh - Calculate genomic callability masks for 1000G and Neanderthal
# ============================================================================

set -e  # Stop immediately if any command fails

JSON=$1
CORES=$2

# Validate input arguments
if [[ -z "$JSON" || -z "$CORES" ]]; then
    echo " Usage: $0 <config.json> <cores>"
    exit 1
fi

# Load parameters and setup paths
PREFIX="$(jq -r '.prefix' "${JSON}")"
mkdir -p "$PREFIX"

# Handle paths with tilde expansion
SIZES_FILE_RAW="$(jq -r '.files.chr_lengths' "$JSON")"
SIZES_FILE="${SIZES_FILE_RAW/#\~/$HOME}"
STRICT_MASK_RAW="$(jq -r '.files["1000GP_files"].bed' "$JSON")"
STRICT_MASK_PATH="${STRICT_MASK_RAW/#\~/$HOME}"

# Get output filenames from JSON
FILE_COV_1KG_NAME="$(jq -r '.window_callability.Thousand_genomes' "${JSON}")"
FILE_COV_ND_NAME="$(jq -r '.window_callability.Nd_1k_genomes' "${JSON}")"

OUT_COV_1KG="$PREFIX/$FILE_COV_1KG_NAME"
OUT_COV_ND="$PREFIX/$FILE_COV_ND_NAME"
TEMP_MASK="$PREFIX/temp_intersect_mask.bed"  # Temporary file for mask intersection

# Performance settings for sort
SORT_ARGS="-S 4G --parallel=$CORES"
export LC_ALL=C

# --- Chromosome Name Normalization ---
CHROM_RAW=$(jq -r '.CHROM' "$JSON")
if [[ ! -f "$SIZES_FILE" ]]; then
    echo "❌ Error: Sizes file not found at $SIZES_FILE"
    exit 1
fi

# Determine chromosome naming convention
FIRST_CHROM_IN_FILE=$(head -1 "$SIZES_FILE" | cut -f1)

# Map JSON chromosome name to sizes file format
if [[ "$FIRST_CHROM_IN_FILE" == chr* && "$CHROM_RAW" != chr* ]]; then
    CHROM_TARGET="chr$CHROM_RAW"
elif [[ "$FIRST_CHROM_IN_FILE" != chr* && "$CHROM_RAW" == chr* ]]; then
    CHROM_TARGET="${CHROM_RAW#chr}"
else
    CHROM_TARGET="$CHROM_RAW"
fi

echo " Processing Chromosome: $CHROM_RAW (Mapped to '$CHROM_TARGET')"
echo " 1kG Output:  $OUT_COV_1KG"
echo " Neand Output: $OUT_COV_ND"

# Check required tools
for cmd in jq bedtools sort zcat; do
    if ! command -v "$cmd" &> /dev/null; then
        echo " Error: '$cmd' not found."
        exit 1
    fi
done

# --- STEP 1: Create Merged Mask ---
# (Neanderthal_1 ∪ Neanderthal_2) ∩ (1000GP Strict Mask)
echo "  Generating merged mask..."

jq -r '.files.neand_files[].bed' "$JSON" | sed "s|^~|$HOME|" \
| xargs -I {} zcat -f {} \
| grep -v "^#" \
| sort $SORT_ARGS -k1,1 -k2,2n \
| bedtools merge -i - \
| bedtools intersect -sorted \
    -a <(zcat -f "$STRICT_MASK_PATH" | grep -v "^#" | sort $SORT_ARGS -k1,1 -k2,2n) \
    -b - > "$TEMP_MASK"

# --- STEP 2: Filter Chromosome Size ---
echo " Filtering chromosome size..."
awk -v c="$CHROM_TARGET" '$1 == c' "$SIZES_FILE" > temp_current.size

if [ ! -s temp_current.size ]; then
    echo " Error: Chromosome '$CHROM_TARGET' not found in sizes file."
    exit 1
fi

# --- STEP 3: Calculate Callability (Coverage) ---
# Calculate coverage in 1000bp windows

echo " Calculating Neanderthal intersection coverage..."
bedtools makewindows -g temp_current.size -w 1000 \
| bedtools coverage -a stdin -b "$TEMP_MASK" > "$OUT_COV_ND"

echo " Calculating 1kG coverage..."
bedtools makewindows -g temp_current.size -w 1000 \
| bedtools coverage -a stdin -b "$STRICT_MASK_PATH" > "$OUT_COV_1KG"

# --- Cleanup ---
rm -f temp_current.size "$TEMP_MASK"

echo ""
echo " Done. Callability files created:"
echo "   - $OUT_COV_ND"
echo "   - $OUT_COV_1KG"

