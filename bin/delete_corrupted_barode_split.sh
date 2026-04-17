#!/bin/bash
# delete_corrupted_barcode_split.sh
#
# Deletes all Nextflow work directories associated with corrupted
# barcode-split cells. Accepts cell IDs as arguments or from a file
# via --file/-f flag.
#
# Usage:
#   # Single or multiple cell IDs as arguments
#   ./delete_corrupted_barcode_split.sh CELL_ID1 CELL_ID2 ...
#
#   # Cell IDs from a file (one per line)
#   ./delete_corrupted_barcode_split.sh -f corrupted_cells.txt
#
# Example:
#   ./delete_corrupted_barcode_split.sh WD1544P_3A_1048_IGO_10626_M_2_b138
#   ./delete_corrupted_barcode_split.sh -f corrupted_cells.txt

set -euo pipefail

usage() {
    echo "Usage: $(basename "$0") [-f file] [CELL_ID ...]"
    echo "  -f FILE   Read cell IDs from FILE (one per line)"
    echo "  CELL_ID   One or more cell IDs as arguments"
    exit 1
}

delete_cell() {
    local cell="$1"
    echo "=== Deleting work directories for: $cell ==="

    find work/ \( \
        -name "${cell}_R*.fastq.gz" -o \
        -name "${cell}.bam" -o \
        -name "${cell}*.bai" -o \
        -name "${cell}*.check*" \
    \) | xargs -I {} dirname {} | sort -u | \
    while read -r workdir; do
        echo "  Deleting: $workdir"
        rm -rf "$workdir"
    done
}

# Parse arguments
CELLS=()

if [ $# -eq 0 ]; then
    usage
fi

while [ $# -gt 0 ]; do
    case "$1" in
        -f|--file)
            [ -f "$2" ] || { echo "ERROR: File not found: $2"; exit 1; }
            while IFS= read -r line; do
                [ -n "$line" ] && CELLS+=("$line")
            done < "$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            CELLS+=("$1")
            shift
            ;;
    esac
done

# Verify work directory exists
[ -d "work/" ] || { echo "ERROR: work/ directory not found"; exit 1; }

# Process cells
echo "Processing ${#CELLS[@]} cell(s)..."
for cell in "${CELLS[@]}"; do
    delete_cell "$cell"
done

echo ""
echo "Deletion complete. Resume with: nextflow run main.nf -resume"
