// FASTQ Quality Control Processes

process FASTQ_SCREEN {
    tag "${cell_id}"
    // no need to publish
    // publishDir "${params.outdir}/results/qc/fastq_screen", mode: 'copy',
    //           pattern: "*_screen.*"

    input:
    tuple val(cell_id), path(fastq_r1), path(fastq_r2)

    output:
    tuple val(cell_id), path("${cell_id}_R1_screen.txt"), 
        path("${cell_id}_R1_screen.html"), emit: screen_results, optional: true
    path "${cell_id}_R1_screen.png", emit: screen_plot, optional: true

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

// // estimate complexity curve
// process PRESEQ_CCURVE {
//     tag "${meta.id}"
//     publishDir "${params.outdir}/results/qc/preseq/c_curve", mode: 'copy'
// 
//     input:
//     tuple val(meta), path(bam), path(bai)
// 
//     output:
//     tuple val(meta), path("${prefix}_c_curve.txt"), emit: ccurve
//     path "versions.yml", emit: versions
// 
//     script:
//     prefix = task.ext.prefix ?: "${meta.id}"
//     def args = task.ext.args ?: ""
//     
//     """
//     /usersoftware/singers/miniforge3/envs/preseq/bin/preseq c_curve \\    
//         -B \\
//         -P \\
//         -o ${prefix}_c_curve.txt \\
//         ${args} \\
//         ${bam}
// 
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         preseq: \$(/usersoftware/singers/miniforge3/envs/preseq/bin/preseq 2>&1 | grep Version | sed 's/Version: //')
//     END_VERSIONS
//     """
// }
// 
// // estimate genome covearage
// process PRESEQ_GCEXTRAP {
//     tag "${meta.id}"
//     publishDir "${params.outdir}/results/qc/preseq/gc_extrap", mode: 'copy'
// 
//     input:
//     tuple val(meta), path(bam), path(bai)
// 
//     output:
//     tuple val(meta), path("${prefix}_gc_extrap.txt"), emit: gcextrap
//     path "versions.yml", emit: versions
// 
//     script:
//     prefix = task.ext.prefix ?: "${meta.id}"
//     def args = task.ext.args ?: ""
//     
//     """
//     /usersoftware/singers/miniforge3/envs/preseq/bin/preseq gc_extrap \\
//         -o ${prefix}_gc_extrap.txt \\
//         ${args} \\
//         ${bam}
// 
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         preseq: \$(/usersoftware/singers/miniforge3/envs/preseq/bin/preseq 2>&1 | grep Version | sed 's/Version: //')
//     END_VERSIONS
//     """
// }
// 
// 
// // estimate read extrapolation
// process PRESEQ_LCEXTRAP {
//     tag "${meta.id}"
//     publishDir "${params.outdir}/results/qc/preseq/lc_extrap", mode: 'copy'
// 
//     input:
//     tuple val(meta), path(bam), path(bai)
// 
//     output:
//     tuple val(meta), path("${prefix}_lc_extrap.txt"), emit: lcextrap
//     path "versions.yml", emit: versions
// 
//     script:
//     prefix = task.ext.prefix ?: "${meta.id}"
//     def args = task.ext.args ?: ""
//     
//     """
//     /usersoftware/singers/miniforge3/envs/preseq/bin/preseq lc_extrap \\
//         -B \\
//         -P \\
//         -o ${prefix}_lc_extrap.txt \\
//         ${args} \\
//         ${bam}
// 
//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         preseq: \$(/usersoftware/singers/miniforge3/envs/preseq/bin/preseq 2>&1 | grep Version | sed 's/Version: //')
//     END_VERSIONS
//     """
// }

