#!/bin/bash
set -euo pipefail

# Parse arguments
INPUT_DIR=""
RESOLUTION=""
OUTPUT_PREFIX="."

while [[ $# -gt 0 ]]; do
    case $1 in
        --input-dir) INPUT_DIR="$2"; shift 2 ;;
        --resolution) RESOLUTION="$2"; shift 2 ;;
        --output-prefix) OUTPUT_PREFIX="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

echo "Processing ${RESOLUTION}k resolution seg files in ${INPUT_DIR}"

# Find seg files
SEG_FILES=($(find "${INPUT_DIR}" -name "*${RESOLUTION}k.bwa_seg.txt" | sort))
NUM_FILES=${#SEG_FILES[@]}

if [ $NUM_FILES -eq 0 ]; then
    echo "No seg files found"
    exit 1
fi

echo "Found $NUM_FILES seg files"

# Create chromosome info from first file
head -1 "${SEG_FILES[0]}" > "${OUTPUT_PREFIX}/chrom_info_${RESOLUTION}k.txt"
tail -n +2 "${SEG_FILES[0]}" | cut -f1-5 >> "${OUTPUT_PREFIX}/chrom_info_${RESOLUTION}k.txt"

echo "Wrote chrom_info_${RESOLUTION}k.txt"

# Extract cell IDs for headers
CELL_IDS=()
for file in "${SEG_FILES[@]}"; do
    cell_id=$(basename "$file" | sed 's/_seg.txt//' | cut -d'.' -f1)
    CELL_IDS+=("$cell_id")
done

# Create matrices using paste (memory efficiency + gzip compressed output)
declare -a COLUMNS=(4 7 8 9 10 11)  # count ratio lowess.ratio seg.mean lowess.ratio.quantal seg.mean.quantal
declare -a NAMES=("counts" "ratio" "lowess_ratio" "seg_mean" "lowess_ratio_quantal" "seg_mean_quantal")

for i in "${!COLUMNS[@]}"; do
    col=${COLUMNS[$i]}
    name=${NAMES[$i]}
    
    echo "Processing ${name} matrix (column ${col})..."
    
    # Create header and pipe to gzip
    echo ${CELL_IDS[@]} | tr ' ' "\t" | gzip -c > "${OUTPUT_PREFIX}/matrix_${name}_${RESOLUTION}k.txt.gz"
    
    # Use paste with process substitution for memory efficiency, pipe to gzip
    paste_cmd=""
    for file in "${SEG_FILES[@]}"; do
        paste_cmd+=" <(tail -n +2 '$file' | cut -f$col)"
    done
    
    eval "paste $paste_cmd" | gzip -c >> "${OUTPUT_PREFIX}/matrix_${name}_${RESOLUTION}k.txt.gz"
    
    echo "Wrote matrix_${name}_${RESOLUTION}k.txt.gz"
done

echo "Matrix aggregation completed"
