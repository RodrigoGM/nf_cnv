// Process to check FASTQ pair synchronization
process CHECK_FASTQ_PAIRS {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(r1_file), path(r2_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_check.tsv"), emit: check_results
    
    script:
    """
    ## Run check_fastq_pairs tool and capture output
    check_fastq_pairs ${r1_file} ${r2_file} > ${sample_id}_check.tsv
    """
}

// Aggregate FASTQ pair check results into summary table
process AGGREGATE_FASTQ_CHECKS {
    publishDir "${params.outdir}/results/qc", mode: 'copy'
    
    input:
    path check_files
    
    output:
    path "fastq_pairs_summary.tsv"
    
    script:
    """
    cat ${check_files} | sort -k1V > fastq_pairs_summary.tsv
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
    tuple val(sample_id), val(barcode_name), path("${sample_id}_${barcode_name}_stats.tsv"), emit: stats

    script:
    """
    # Execute barcode_split with default parameters
    barcode_split \\
        -p ${barcode_seq} \\
        -u "TGTGTTGGGTGTGTTTGG" \\
        -o ${sample_id}_${barcode_name} \\
        ${r1_file} ${r2_file}

    # Ensure output files exist (create empty gzipped files if needed)
    if [[ ! -f "${sample_id}_${barcode_name}_R1.fastq.gz" ]]; then
        echo | gzip > ${sample_id}_${barcode_name}_R1.fastq.gz
    fi
    if [[ ! -f "${sample_id}_${barcode_name}_R2.fastq.gz" ]]; then
        echo | gzip > ${sample_id}_${barcode_name}_R2.fastq.gz
    fi

    # Count reads - check if files have content
    if [[ -s "${sample_id}_${barcode_name}_R1.fastq.gz" ]]; then
        read_count=\$(zcat ${sample_id}_${barcode_name}_R1.fastq.gz | wc -l | awk '{print \$1/4}')
    else
        read_count=0
    fi

    file_size=\$(du -h ${sample_id}_${barcode_name}_R1.fastq.gz | cut -f1)
    md5=\$( md5sum ${sample_id}_${barcode_name}_R1.fastq.gz )

    echo -e "\$md5\\t${sample_id}\\t${barcode_name}\\t${barcode_seq}\\t\$read_count\\t\$file_size" > ${sample_id}_${barcode_name}_stats.tsv
    """
}

// check fastq pairs
process CHECK_CELL_FASTQ_PAIRS {
    tag "${cell_id}"
    
    input:
    tuple val(cell_id), path(r1_file), path(r2_file)
    
    output:
    path "${cell_id}_check.tsv"
    
    script:
    """
    check_fastq_pairs ${r1_file} ${r2_file} > ${cell_id}_check.tsv
    """
}

// Aggregate cell-level FASTQ check results
process AGGREGATE_CELL_CHECKS {
    publishDir "${params.outdir}/results/qc", mode: 'copy'
    
    input:
    path check_files
    
    output:
    path "cell_fastq_pairs_summary.tsv"
    
    script:
    """
    cat ${check_files} | sort -k1V > cell_fastq_pairs_summary.tsv
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
    ## Create TSV with header
    echo -e "md5\\tfile_name\\tsample_id\\tbarcode_name\\tbarcode_seq\\tcell_reads\tfile_size" > barcode_split_summary.tsv
    
    ## Concatenate all statistics files
    cat ${stats_files} | sort -k1V -k2V >> barcode_split_summary.tsv
    """
}
