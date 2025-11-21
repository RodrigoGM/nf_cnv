// CNV Analysis using varbin from cna_utils
process GET_BIN_COUNTS {
    tag "${cell_id}_${resolution}k"
    // conditional publishDir in modules.config
    // publishDir "${params.outdir}/results/varbin${resolution}k/counts", mode: 'copy',
    //            pattern: "*.bin.counts.bed"
    // publishDir "${params.outdir}/results/varbin${resolution}k/stats", mode: 'copy',
    //            pattern: "*.bin.counts.stats.bed"

    input:
    tuple val(cell_id), path(strand_bam), path(strand_bai), val(resolution)

    output:
    tuple val(cell_id), val(resolution), 
          path("${cell_id}.${params.genome}.PE_*.${resolution}k.bwa.bin.counts.bed"), 
          path("${cell_id}.${params.genome}.PE_*.${resolution}k.bwa.bin.counts.stats.bed"), emit: bin_counts

    script:
    def genome_map = [hsa37: 'hg19', hsa38: 'hg38']
    def cna_genome = genome_map[params.genome]
    def strand_ext = params.use_reverse_reads ? 'PE_RV' : 'PE_FW'
    def output_prefix = "${cell_id}.${params.genome}.${strand_ext}.${resolution}k.bwa"
    
    """
    echo "${resolution}k bin count for ${cell_id}"
    
    ${params.cna_utils_path}/scripts/getBinCounts.py \\
        -i ${strand_bam} \\
        -b ${params.cna_utils_path}/data/${cna_genome}_${resolution}k_gz_enc_bins.bed \\
        -d ${params.cna_utils_path}/data/${cna_genome}_150bp_dz_enc.bed \\
        -o ${output_prefix}.bin.counts.bed \\
        -v > ${output_prefix}.bin.counts.stats.bed
    """
}

// create different directories for each output
process CNV_PROFILE {
    tag "${cell_id}_${resolution}k"
    
    // publishDir "${params.outdir}/results/varbin${resolution}k/seg_long", mode: 'copy',
    //           pattern: "*_seg.txt"
    // publishDir "${params.outdir}/results/varbin${resolution}k/seg_short", mode: 'copy',
    //           pattern: "*_short_seg.txt"
    // publishDir "${params.outdir}/results/varbin${resolution}k/ploidy", mode: 'copy',
    //           pattern: "*.quantal.ploidy.txt"
    // publishDir "${params.outdir}/results/varbin${resolution}k/logs", mode: 'copy',
    //           pattern: "*.quantal.log"

    input:
    tuple val(cell_id), val(resolution), path(bin_counts), path(bin_stats)

    output:
    tuple val(cell_id), val(resolution),
          path("${cell_id}.${params.genome}.PE_*.${resolution}k.bwa.quantal.ploidy.txt"), emit: ploidy_results
    path "${cell_id}.${params.genome}.PE_*.${resolution}k.bwa.quantal.log", emit: log_files
    path "${cell_id}.${params.genome}.PE_*.${resolution}k.bwa_seg.txt", emit: seg_files
    path "${cell_id}.${params.genome}.PE_*.${resolution}k.bwa_short_seg.txt", emit: short_seg_files
    path "${cell_id}.${params.genome}.PE_*.${resolution}k.bwa*", emit: all_results

    script:
    def genome_map = [hsa37: 'hg19', hsa38: 'hg38']
    def cna_genome = genome_map[params.genome]
    def strand_ext = params.use_reverse_reads ? 'PE_RV' : 'PE_FW'
    def output_prefix = "${cell_id}.${params.genome}.${strand_ext}.${resolution}k.bwa"

    """
    echo "${resolution}k varbin CN estimation for ${cell_id}"

    ${params.cna_utils_path}/scripts/cnvProfile.R \\
        -b ${bin_counts} \\
        -g ${params.cna_utils_path}/data/${cna_genome}_${resolution}k_gz_enc_gc.bed \\
        -n ${output_prefix} \\
        --minploidy=${params.min_ploidy} \\
        --maxploidy=${params.max_ploidy} \\
        -v 2> ${output_prefix}.quantal.log

    # Extract ploidy and error information
    echo -e "cellID\\tploidy\\terror" > ${output_prefix}.quantal.ploidy.txt

    ploidy=\$(grep "Ploidy" ${output_prefix}.quantal.log | tail -n 1 | sed -e 's/.*Ploidy: //' | awk '{printf "%.4f", \$1}' || echo "NA")
    error=\$(grep "Error" ${output_prefix}.quantal.log | tail -n 1 | sed -e 's/.*Error: //' | awk '{printf "%.4f", \$1}' || echo "NA")

    echo -e "${cell_id}\\t\${ploidy}\\t\${error}" >> ${output_prefix}.quantal.ploidy.txt
    
    # Create seg files if they don't exist
    for suffix in _seg.txt _short_seg.txt; do
        if [[ ! -f "${output_prefix}\${suffix}" ]]; then
            touch "${output_prefix}\${suffix}"
        fi
    done
    """
}

// Create ploidy files
process AGGREGATE_CNV_PLOIDY {
    tag "all_samples"
    publishDir "${params.outdir}/results/cnv_summary", mode: 'copy'

    input:
    path ploidy_files

    output:
    path "ploidy_results_combined.tsv", emit: combined_results
    path "ploidy_results_*.tsv", emit: resolution_results

    script:
    """
    # Create combined results file
    echo -e "cellID\\tploidy\\terror\\tresolution\\tgenome\\tread" > ploidy_results_combined.tmp
    
    # Process each ploidy file
    for file in ${ploidy_files}; do
        # Extract resolution from filename (format: ...5k.bwa.quantal.ploidy.txt)
        if [[ \$file == *".5k.bwa.quantal.ploidy.txt" ]]; then
            resolution="5k"
        elif [[ \$file == *".20k.bwa.quantal.ploidy.txt" ]]; then
            resolution="20k"
        elif [[ \$file == *".50k.bwa.quantal.ploidy.txt" ]]; then
            resolution="50k"
        else
            resolution="unknown"
        fi
        
        # Extract read type (FW/RV) from filename
        if [[ \$file == *"PE_FW"* ]]; then
            read_type="FW"
        elif [[ \$file == *"PE_RV"* ]]; then
            read_type="RV"
        else
            read_type="unknown"
        fi
        
        # Skip header and add resolution/genome/read info
        if [[ -f \$file ]]; then
            tail -n +2 \$file | while IFS='\\t' read -r cellID ploidy error; do
                if [[ -n "\$cellID" ]]; then
                    echo -e "\${cellID}\\t\${ploidy}\\t\${error}\\t\${resolution}\\t${params.genome}\\t\${read_type}" >> ploidy_results_combined.tmp
                fi
            done
        fi
    done
    
    cat ploidy_results_combined.tmp | tr "\t" ' ' | sed -E 's/ +/ /g' | tr ' ' "\t" > ploidy_results_combined.tsv
    
    """
}

process CNV_PLOIDY_SUMMARY {
    tag "summary_stats"
    publishDir "${params.outdir}/results/cnv_summary", mode: 'copy'

    input:
    path combined_results

    output:
    path "ploidy_summary_all.txt", emit: ploidy_stats

    script:
    """
    # Run summary statistics with grouping by resolution
    ${projectDir}/bin/cnv_ploidy_summary.py \\
        -i ${combined_results} \\
        -o ploidy_summary_all.txt \\
        --group-by resolution
  
    """
}

process AGGREGATE_SEG_MATRICES {
    tag "all_samples_${resolution}k"
    publishDir "${params.outdir}/results/cnv_summary/", mode: 'copy'

    input:
    tuple val(resolution), path(seg_files)

    output:
    path "chrom_info_${resolution}k.txt", emit: chrom_info
    path "matrix_counts_${resolution}k.txt.gz", emit: counts_matrix
    path "matrix_ratio_${resolution}k.txt.gz", emit: ratio_matrix
    path "matrix_lowess_ratio_${resolution}k.txt.gz", emit: lowess_ratio_matrix
    path "matrix_seg_mean_${resolution}k.txt.gz", emit: seg_mean_matrix
    path "matrix_lowess_ratio_quantal_${resolution}k.txt.gz", emit: lowess_ratio_quantal_matrix
    path "matrix_seg_mean_quantal_${resolution}k.txt.gz", emit: seg_mean_quantal_matrix

    script:
    """
    ${projectDir}/bin/aggregate_seg_matrices.sh \\
        --input-dir . \\
        --resolution ${resolution} \\
        --output-prefix .
    """
}

process AGGREGATE_BIN_COUNTS_STATS {
    tag "bin_counts_stats_${resolution}k"
    publishDir "${params.outdir}/results/cnv_summary", mode: 'copy'

    input:
    tuple val(resolution), path(stats_files)

    output:
    path "bin_counts_stats_${resolution}k.txt", emit: agg_stats

    script:
    """
    echo -e "cellID\\ttotal.aligned\\tpassed.reads\\tfiltered.reads\\tcounted.reads" > bin_counts_stats_${resolution}k.txt
    for file in ${stats_files}; do
        cell_id=\$(basename \$file | cut -d. -f1)
        total=\$(awk '/Total alignments:/ {print \$3}' \$file)
        passed=\$(awk '/Passed alignments:/ {print \$3}' \$file)
        filtered=\$(awk '/Filtered alignments:/ {print \$3}' \$file)
        counted=\$(awk '/Counted alignments:/ {print \$3}' \$file)
        echo -e "\${cell_id}\\t\${total}\\t\${passed}\\t\${filtered}\\t\${counted}" >> bin_counts_stats_${resolution}k.txt
    done
    """
}

process CREATE_CELL_PHENOTYPE_TEMPLATE {
    tag "phenotype_template"
    publishDir "${params.outdir}/results/cnv_summary", mode: 'copy'

    input:
    path ploidy_20k_files

    output:
    path "cell_phenotype.txt", emit: phenotype_template

    script:
    """
    # Create header
    echo -e "cellID\\tbioID\\tsampleID\\tsubsample\\tplate\\tsubtype\\tmold.histology\\tmold.grade\\tmold.viability\\tmold.appearance\\tbarcode\\tsample.type\\tgate\\tpct.[HISTOLOGY1]\\tpct.[HISTOLOGY...]\\tpath.comments\\tploidy\\terror" > cell_phenotype.txt
    
    # Extract cell data from 20k ploidy files
    for file in ${ploidy_20k_files}; do
        if [[ \$file == *"20k"* ]]; then
            tail -n +2 \$file | while IFS='\t' read -r cellID ploidy error; do
                if [[ -n "\$cellID" ]]; then
                    echo -e "\${cellID}\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\\t\${ploidy}\\t\${error}" >> cell_phenotype.txt
                fi
            done
        fi
    done
    
    echo "Created phenotype template with \$(tail -n +2 cell_phenotype.txt | wc -l) cells"
    """
}
