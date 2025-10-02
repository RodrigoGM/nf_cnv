// Simple file discovery process
process DISCOVER_FASTQ_FILES {
    publishDir "${params.outdir}/results/merged", mode: 'copy'
    
    output:
    path "fastq_manifest.txt"
    stdout
    
    script:
    """
    if [[ "${params.sample_pattern}" != ".*" ]]; then
        R1=( \$(find ${params.seqdata} -name "*${params.sample_pattern}*_R1_001.fastq.gz") )
    else
        R1=( \$(find ${params.seqdata} -name "*_R1_001.fastq.gz") )
    fi
    
    R2=( \$(for i in "\${R1[@]}"; do echo "\$i" | sed -e 's/_R1_001/_R2_001/'; done) )
    sample_id=( \$(for i in "\${R1[@]}"; do echo "\$i" | sed -e 's/_IGO.*//'; done) )
    
    nn=\$(( \${#R1[@]} - 1 ))
    for i in \$(seq 0 \$nn); do
        echo -e "\${sample_id[i]}\t\${R1[i]}\t\${R2[i]}"
    done > fastq_manifest.txt
    
    # Output for channel processing
    for i in \$(seq 0 \$nn); do
        echo "\${sample_id[i]}|\${R1[i]}|\${R2[i]}"
    done
    """
}


// Process to organize files by sample for CAT_FASTQ
process ORGANIZE_FOR_CONCAT {
    tag "organize"
    
    input:
    path manifest
    
    output:
    path "sample_*.txt", emit: sample_lists
    
    script:
    """
    #!/bin/bash
    
    # Get unique sample IDs
    tail -n +2 ${manifest} | cut -d',' -f1 | sort | uniq > samples.txt
    
    # Create file lists for each sample
    while read -r sample_id; do
        echo "Processing sample: \$sample_id"
        
        # Get R1 files for this sample
        grep "^\$sample_id,R1," ${manifest} | cut -d',' -f3 | sort -V > "sample_\${sample_id}_R1.txt"
        
        # Get R2 files for this sample  
        grep "^\$sample_id,R2," ${manifest} | cut -d',' -f3 | sort -V > "sample_\${sample_id}_R2.txt"
        
        # Combine for CAT_FASTQ (R1 files first, then R2 files)
        cat "sample_\${sample_id}_R1.txt" "sample_\${sample_id}_R2.txt" > "sample_\${sample_id}.txt"
        
        echo "Sample \$sample_id: \$(wc -l < sample_\${sample_id}_R1.txt) R1 files, \$(wc -l < sample_\${sample_id}_R2.txt) R2 files"
        
    done < samples.txt
    
    echo "Created file lists for \$(wc -l < samples.txt) samples"
    ls -la sample_*.txt
    """
}
