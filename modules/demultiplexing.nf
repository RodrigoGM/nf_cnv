process CHECK_FASTQ_PAIRS {
    tag "${sample_id}"
    
    input:
    tuple val(sample_id), path(r1_file), path(r2_file)
    
    output:
    tuple val(sample_id), path("${sample_id}_check.txt"), emit: check_results
    
    script:
    """
    /usersoftware/singers/opt/bin/check_fastq_pairs ${r1_file} ${r2_file} > ${sample_id}_check.txt
    """
}

process AGGREGATE_FASTQ_CHECKS {
    publishDir "${params.outdir}/results/qc", mode: 'copy'
    
    input:
    path check_files
    
    output:
    path "fastq_pairs_summary.tsv"
    
    script:
    """
    echo -e "sample_id\\tr1_file\\tr1_reads\\tr2_reads\\tsynchronized" > fastq_pairs_summary.tsv
    
    for check_file in ${check_files}; do
        sample=\$(basename "\$check_file" _check.txt)
        while read -r line; do
            echo -e "\$sample\\t\$line" >> fastq_pairs_summary.tsv
        done < "\$check_file"
    done
    """
}
