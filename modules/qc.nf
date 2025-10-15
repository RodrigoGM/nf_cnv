// FASTQ Quality Control Processes

process FASTQ_SCREEN {
    tag "${cell_id}"
    publishDir "${params.outdir}/results/qc/fastq_screen", mode: 'copy',
               pattern: "*_screen.*"

    input:
    tuple val(cell_id), path(fastq_r1), path(fastq_r2)

    output:
    tuple val(cell_id), path("${cell_id}_R1_screen.txt"), 
          path("${cell_id}_R1_screen.html"), emit: screen_results
    path "${cell_id}_R1_screen.png", optional: true, emit: screen_plot

    script:
    def subset_flag = params.fastq_screen_subset > 0 ? 
        "--subset ${params.fastq_screen_subset}" : ""
    
    """
    echo "Running FastQ Screen for ${cell_id} (R1 only)"
    
    fastq_screen \\
        --conf ${params.fastq_screen_config} \\
        --threads ${task.cpus} \\
        ${subset_flag} \\
        --outdir . \\
        ${fastq_r1}
    
    # Rename output files to include cell_id
    if [ -f "*_screen.txt" ]; then
        mv *_screen.txt ${cell_id}_R1_screen.txt
    fi
    if [ -f "*_screen.html" ]; then
        mv *_screen.html ${cell_id}_R1_screen.html
    fi
    if [ -f "*_screen.png" ]; then
        mv *_screen.png ${cell_id}_R1_screen.png
    fi
    """
}

