#!/bin/bash
# delete_corrupted_barcode_split.sh

# For each corrupted cell, find and delete ALL associated work directories
while read cell; do
    echo "=== Deleting work directories for: $cell ==="
    
    # Find all work directories containing this cell's files
    find work/ -name "${cell}_R*.fastq.gz" -o -name "${cell}.bam" -o -name "${cell}*.bai" -o -name "${cell}*.check*" | \
        xargs -I {} dirname {} | \
        sort -u | \
        while read workdir; do
            echo "  Deleting: $workdir"
            rm -rf "$workdir"
        done
    
done < $1

echo ""
echo "Deletion complete. Now run:"
echo "nextflow run main.nf -resume"
