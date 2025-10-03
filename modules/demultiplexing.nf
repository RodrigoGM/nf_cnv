// Process to check FASTQ pair synchronization
process CHECK_FASTQ_PAIRS {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(r1_file), path(r2_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_check.txt"), emit: check_results
    
    script:
    """
    // Run check_fastq_pairs tool and capture output
    check_fastq_pairs ${r1_file} ${r2_file} > ${sample_id}_check.txt
    """
}

// Aggregate FASTQ pair check results into summary table
process AGGREGATE_FASTQ_CHECKS {
    publishDir "${params.outdir}/results/qc", mode: 'copy'
    
    input:
    path check_files
    
    output:
    path "fastq_pairs_summary.txt"
    
    script:
    """
    cat ${check_files} | sort -k1V > cell_fastq_pairs_summary.txt
    """
}

// Split merged FASTQ files by barcode sequences
process BARCODE_SPLIT {
    tag "${sample_id}_${barcode_name}"
    publishDir "${params.outdir}/results/wsplit", mode: 'copy', pattern: "*_R{1,2}.fastq.gz"
    
    input:
    tuple val(sample_id), path(r1_file), path(r2_file), val(barcode_name), val(barcode_seq)
    
    output:
    tuple val("${sample_id}_${barcode_name}"), path("${sample_id}_${barcode_name}_R1.fastq.gz"), path("${sample_id}_${barcode_name}_R2.fastq.gz"), emit: cell_fastqs
    tuple val(sample_id), val(barcode_name), path("${sample_id}_${barcode_name}_stats.txt"), emit: stats
    
    script:
    """
    // Execute barcode_split with standard parameters
    barcode_split \\
        -p ${barcode_seq} \\
        -m 1 \\
        -l 80 \\
        -o ${sample_id}_${barcode_name} \\
        ${r1_file} ${r2_file}
    
    // Count reads in output files using check_fastq_pairs
    if [[ -f "${sample_id}_${barcode_name}_R1.fastq.gz" && -f "${sample_id}_${barcode_name}_R2.fastq.gz" ]]; then
        read_count=\$(check_fastq_pairs ${sample_id}_${barcode_name}_R1.fastq.gz ${sample_id}_${barcode_name}_R2.fastq.gz | cut -f2)
    else
        read_count=0
    fi
    
    // Create statistics file
    echo -e "${sample_id}\\t${barcode_name}\\t${barcode_seq}\\t\$read_count" > ${sample_id}_${barcode_name}_stats.txt
    """
}

// Check cell-level FASTQ pair synchronization
process CHECK_CELL_FASTQ_PAIRS {
    tag "${cell_id}"
    
    input:
    tuple val(cell_id), path(r1_file), path(r2_file)
    
    output:
    path "${cell_id}_check.txt"
    
    script:
    """
    check_fastq_pairs ${r1_file} ${r2_file} > ${cell_id}_check.txt
    """
}

// Aggregate cell-level FASTQ check results
process AGGREGATE_CELL_CHECKS {
    publishDir "${params.outdir}/results/qc", mode: 'copy'
    
    input:
    path check_files
    
    output:
    path "cell_fastq_pairs_summary.txt"
    
    script:
    """
    cat ${check_files} | sort -k1V > cell_fastq_pairs_summary.txt
    """
}


// Aggregate barcode splitting statistics
process AGGREGATE_BARCODE_STATS {
    publishDir "${params.outdir}/results/qc", mode: 'copy'
    
    input:
    path stats_files
    
    output:
    path "barcode_split_summary.tsv"
    
    script:
    """
    // Create TSV with header
    echo -e "sample_id\\tbarcode_name\\tbarcode_seq\\tcell_reads" > barcode_split_summary.tsv
    
    // Concatenate all statistics files
    cat ${stats_files} | sort -k1V -k2V >> barcode_split_summary.tsv
    """
}
