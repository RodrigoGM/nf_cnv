// Initial BWA-MEM alignment 
process BWA_MEM {
    tag "${cell_id}"
    
    input:
    tuple val(cell_id), path(r1_file), path(r2_file)
    
    output:
    tuple val(cell_id), path("${cell_id}.unsorted.bam"), emit: unsorted_bam
    
    script:
    def bwa_index = params.genome_indices[params.genome]
    def rg_line = "@RG\\tID:${cell_id}\\tSM:${cell_id}\\tPL:ILLUMINA\\tLB:${cell_id}"
    """
    # Check if FASTQ files have content
    r1_reads=\$(zcat ${r1_file} | head -n 4 | wc -l)
    r2_reads=\$(zcat ${r2_file} | head -n 4 | wc -l)
    
    if [[ \$r1_reads -eq 0 || \$r2_reads -eq 0 ]]; then
        echo "Empty FASTQ files detected, creating empty BAM"
        bwa mem -aM -t ${task.cpus} -R "${rg_line}" \\
            ${bwa_index} <(echo | gzip -c) <(echo | gzip -c) | \\
            samtools view -bS -F 4 - > ${cell_id}.unsorted.bam
    else
        bwa mem -aM -t ${task.cpus} -R "${rg_line}" \\
            ${bwa_index} ${r1_file} ${r2_file} | \\
            samtools view -bS -F 4 - > ${cell_id}.unsorted.bam
    fi
    """
}

// samtools collate / group sorting to fix mate
process COLLATE_BAM {
    tag "${cell_id}"
    
    input:
    tuple val(cell_id), path(unsorted_bam)
    
    output:
    tuple val(cell_id), path("${cell_id}.collated.bam"), emit: collated_bam
    
    script:
    """
    samtools collate -o ${cell_id}.collated.bam ${unsorted_bam}
    """
}

// samtools fix mate information
process FIXMATE_BAM {
    tag "${cell_id}"
    
    input:
    tuple val(cell_id), path(collated_bam)
    
    output:
    tuple val(cell_id), path("${cell_id}.fixmate.bam"), emit: fixmate_bam
    
    script:
    """
    samtools fixmate -m ${collated_bam} ${cell_id}.fixmate.bam
    """
}

// samtools sort bam
process SORT_BAM {
    tag "${cell_id}"
    
    input:
    tuple val(cell_id), path(fixmate_bam)
    
    output:
    tuple val(cell_id), path("${cell_id}.sorted.bam"), emit: sorted_bam
    
    script:
    """
    samtools sort -@ ${task.cpus} -o ${cell_id}.sorted.bam ${fixmate_bam}
    """
}

// Picard mark duplicates without deleting
process MARK_DUPLICATES {
    tag "${cell_id}"
    publishDir "${params.outdir}/results/bwa_out_md", mode: 'copy',
              pattern: "*.${params.genome}.PE.md.bam*"

    input:
    tuple val(cell_id), path(sorted_bam)

    output:
    tuple val(cell_id), path("${cell_id}.${params.genome}.PE.md.bam"),
          path("${cell_id}.${params.genome}.PE.md.bam.bai"), emit: marked_bam
    path "${cell_id}.dup_metrics.txt", emit: markdup_metrics

    script:
    def optical_params = params.no_optical_duplicates ? 
        "-READ_NAME_REGEX null" : 
        "-OPTICAL_DUPLICATE_PIXEL_DISTANCE ${params.optical_duplicate_pixel_distance}"
    
    """
    # Check if BAM has alignments
    read_count=\$(samtools view -c ${sorted_bam})
    
    if [[ \$read_count -eq 0 ]]; then
        echo "No aligned reads, creating empty output"
        cp ${sorted_bam} ${cell_id}.${params.genome}.PE.md.bam
        samtools index ${cell_id}.${params.genome}.PE.md.bam
        echo -e "## No reads for duplicate marking\\nLIBRARY\\tUNPAIRED_READS_EXAMINED\\tREAD_PAIRS_EXAMINED\\tSECONDARY_OR_SUPPLEMENTARY_RDS\\tUNMAPPED_READS\\tUNPAIRED_READ_DUPLICATES\\tREAD_PAIR_DUPLICATES\\tREAD_PAIR_OPTICAL_DUPLICATES\\tPERCENT_DUPLICATION\\tESTIMATED_LIBRARY_SIZE\\n${cell_id}\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0" > ${cell_id}.dup_metrics.txt
    else
        picard MarkDuplicates \\
            -INPUT ${sorted_bam} \\
            -OUTPUT ${cell_id}.${params.genome}.PE.md.bam \\
            -METRICS_FILE ${cell_id}.dup_metrics.txt \\
            -CREATE_INDEX true \\
            -VALIDATION_STRINGENCY LENIENT \\
            -ASSUME_SORT_ORDER coordinate \\
            ${optical_params}

    mv ${cell_id}.${params.genome}.PE.md.bai  ${cell_id}.${params.genome}.PE.md.bam.bai
    fi
    """
}

// Picard remove duplicates
process REMOVE_DUPLICATES {
    tag "${cell_id}"
    publishDir "${params.outdir}/results/bwa_out_dd", mode: 'copy',
              pattern: "*.${params.genome}.PE.dd.bam*"

    input:
    tuple val(cell_id), path(marked_bam), path(marked_bai)

    output:
    tuple val(cell_id), path("${cell_id}.${params.genome}.PE.dd.bam"),
          path("${cell_id}.${params.genome}.PE.dd.bam.bai"), emit: dedup_bam
    path "${cell_id}.dd_metrics.txt", emit: dedup_metrics

    script:
    def optical_params = params.no_optical_duplicates ?
        "-READ_NAME_REGEX null" :
        "-OPTICAL_DUPLICATE_PIXEL_DISTANCE ${params.optical_duplicate_pixel_distance ?: 2500}"
    
    """
    # Check if BAM has alignments
    read_count=\$(samtools view -c ${marked_bam})
    
    if [[ \$read_count -eq 0 ]]; then
        echo "No aligned reads, creating empty output"
        cp ${marked_bam} ${cell_id}.${params.genome}.PE.dd.bam
        samtools index ${cell_id}.${params.genome}.PE.dd.bam
        echo -e "## No reads for duplicate removal\\nLIBRARY\\tUNPAIRED_READS_EXAMINED\\tREAD_PAIRS_EXAMINED\\tSECONDARY_OR_SUPPLEMENTARY_RDS\\tUNMAPPED_READS\\tUNPAIRED_READ_DUPLICATES\\tREAD_PAIR_DUPLICATES\\tREAD_PAIR_OPTICAL_DUPLICATES\\tPERCENT_DUPLICATION\\tESTIMATED_LIBRARY_SIZE\\n${cell_id}\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0\\t0" > ${cell_id}.dd_metrics.txt
    else
        picard MarkDuplicates \\
            -INPUT ${marked_bam} \\
            -OUTPUT ${cell_id}.${params.genome}.PE.dd.bam \\
            -METRICS_FILE ${cell_id}.dd_metrics.txt \\
            -REMOVE_DUPLICATES true \\
            -CREATE_INDEX true \\
            -VALIDATION_STRINGENCY LENIENT \\
            -ASSUME_SORT_ORDER coordinate \\
            ${optical_params}
    mv ${cell_id}.${params.genome}.PE.dd.bai ${cell_id}.${params.genome}.PE.dd.bam.bai
    fi
    """
}

// Quality Filter, and Mate Clipping
process FILTER_AND_EXTRACT_READS {
    tag "${cell_id}"
    publishDir "${params.outdir}/results/${params.use_reverse_reads ? 'bwa_out_rv' : 'bwa_out_fw'}",
	       mode: 'copy',
               pattern: "*.${params.genome}.PE_*.dd.bam*"

    input:
    tuple val(cell_id), path(dedup_bam), path(dedup_bai)

    output:
    tuple val(cell_id), path("${cell_id}.${params.genome}.PE_*.dd.bam"),
          path("${cell_id}.${params.genome}.PE_*.dd.bam.bai"), emit: strand_bam

    script:
    // Define variables at script level so they're available for output declarations
    def strand_flag = params.use_reverse_reads ? '0x80' : '0x40'
    def strand_suffix = params.use_reverse_reads ? 'RV' : 'FW'
    
    """
    # Filter unique reads and extract strand in one step
    samtools view -@ ${task.cpus} -q 30 -F 0x800 -f ${strand_flag} -h -b \\
        -o ${cell_id}.${params.genome}.PE_${strand_suffix}.dd.bam \\
        ${dedup_bam} ${params.primary_chromosomes}
    
    samtools index ${cell_id}.${params.genome}.PE_${strand_suffix}.dd.bam
    """
}

// Utilities
process VALIDATE_BAM {
    tag "${cell_id}"
    
    input:
    tuple val(cell_id), path(bam_file)
    
    output:
    tuple val(cell_id), path(bam_file), emit: validated_bam
    
    script:
    """
    # Quick validation
    samtools quickcheck ${bam_file}
    
    # Check for duplicate read names
    dup_count=\$(samtools view ${bam_file} | cut -f1 | sort | uniq -d | wc -l)
    if [[ \$dup_count -gt 0 ]]; then
        echo "WARNING: Found \$dup_count duplicate read names in ${cell_id}"
    fi
    
    echo "BAM validation passed for ${cell_id}"
    """
}

