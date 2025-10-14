// CNV Analysis using varbin from cna_utils
process GET_BIN_COUNTS {
    tag "${cell_id}_${resolution}k"
    publishDir "${params.outdir}/results/varbin${resolution}k", mode: 'copy',
               pattern: "*.bed"

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

process CNV_PROFILE {
    tag "${cell_id}_${resolution}k"
    publishDir "${params.outdir}/results/varbin${resolution}k", mode: 'copy',
               pattern: "${cell_id}.${params.genome}.${strand_ext}.${resolution}k.bwa.*"

    input:
    tuple val(cell_id), val(resolution), path(bin_counts), path(bin_stats)

    output:
    tuple val(cell_id), val(resolution), 
          path("${cell_id}.${params.genome}.${strand_ext}.${resolution}k.bwa.quantal.ploidy.txt"), emit: ploidy_results
    path "${cell_id}.${params.genome}.${strand_ext}.${resolution}k.bwa.quantal.log", emit: log_files
    path "${cell_id}.${params.genome}.${strand_ext}.${resolution}k.bwa_seg.txt", emit: seg_files
    path "${cell_id}.${params.genome}.${strand_ext}.${resolution}k.bwa_short_seg.txt", emit: short_seg_files
    path "${cell_id}.${params.genome}.${strand_ext}.${resolution}k.bwa*", emit: all_results

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
    
    ploidy=\$(grep "Ploidy" ${output_prefix}.quantal.log | tail -n 1 | sed -e 's/.*Ploidy: //' || echo "NA")
    error=\$(grep "Error" ${output_prefix}.quantal.log | tail -n 1 | sed -e 's/.*Error: //' | awk '{printf "%.4f", \$1}' || echo "NA")
    
    echo -e "${cell_id}\\t\${ploidy}\\t\${error}" >> ${output_prefix}.quantal.ploidy.txt
    """
}

process AGGREGATE_CNV_RESULTS {
    tag "all_samples"
    publishDir "${params.outdir}/results/cnv_summary", mode: 'copy'

    input:
    path ploidy_files

    output:
    path "combined_ploidy_results.txt", emit: combined_results
    path "*_resolution_summary.txt", emit: resolution_summaries

    script:
    """
    # Create combined results file
    echo -e "cellID\\tresolution\\tploidy\\terror\\tgenome\\tstrand" > combined_ploidy_results.txt
    
    # Process each ploidy file
    for file in ${ploidy_files}; do
        if [[ \$file == *"5k"* ]]; then
            resolution="5k"
        elif [[ \$file == *"20k"* ]]; then
            resolution="20k"
        elif [[ \$file == *"50k"* ]]; then
            resolution="50k"
        else
            resolution="unknown"
        fi
        
        strand="${params.use_reverse_reads ? 'RV' : 'FW'}"
        
        # Skip header and add resolution info
        tail -n +2 \$file | while read line; do
            echo -e "\${line}\\t\${resolution}\\t${params.genome}\\t\${strand}" >> combined_ploidy_results.txt
        done
    done
    
    # Create resolution-specific summaries
    for res in 5k 20k 50k; do
        echo "Summary for \${res} resolution:" > \${res}_resolution_summary.txt
        echo "=========================" >> \${res}_resolution_summary.txt
        
        grep "\t\${res}\t" combined_ploidy_results.txt | \\
        awk -F'\\t' 'NR>1 {
            if(\$3 != "NA") {
                sum += \$3; count++; 
                if(\$3 < min || min=="") min = \$3;
                if(\$3 > max || max=="") max = \$3;
            }
        } 
        END {
            if(count > 0) {
                print "Cells processed: " count
                print "Mean ploidy: " sum/count
                print "Min ploidy: " min
                print "Max ploidy: " max
            } else {
                print "No valid ploidy estimates found"
            }
        }' >> \${res}_resolution_summary.txt
    done
    """
}
