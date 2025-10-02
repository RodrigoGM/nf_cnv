#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { CAT_FASTQ } from './modules/nf-core/cat/fastq/main'
include { discoverFastqFiles; createManifestData; discoverMergedFastqPairs } from './modules/utils'
include { CHECK_FASTQ_PAIRS; AGGREGATE_FASTQ_CHECKS } from './modules/demultiplexing'


// Parameters
params.seqdata = "${launchDir}/seqdata"
params.outdir = "${launchDir}"
params.sample_pattern = ".*"

process CREATE_FASTQ_MANIFEST {
    publishDir "${params.outdir}/results/qc", mode: 'copy'
    
    input:
    val manifest_data
    
    output:
    path "fastq_manifest.csv"
    
    script:
    """
    cat > fastq_manifest.csv << 'EOF'
${manifest_data}
EOF
    """
}



workflow {
    log.info "SC-CNV Processing Pipeline"
    
    // Discover FASTQ files by sample
    samples = discoverFastqFiles(params.seqdata, params.sample_pattern)
    log.info "Found ${samples.size()} samples matching pattern '${params.sample_pattern}'"
    
    if (samples.isEmpty()) {
        error "No samples found matching pattern '${params.sample_pattern}' in ${params.seqdata}"
    }
    
    // Create channel for CAT_FASTQ
    sample_ch = Channel.fromList(
        samples.collect { sample_id, files ->
            def meta = [id: sample_id, single_end: false]
            def reads = [files.r1, files.r2].flatten().sort { it.name }
            [meta, reads]
        }
    )
    
    // Concatenate FASTQ files
    CAT_FASTQ(sample_ch)
    
    // Create manifest
    manifest_data = createManifestData(samples)
    CREATE_FASTQ_MANIFEST(Channel.value(manifest_data))
    
    log.info "Concatenation completed for ${samples.size()} samples"


    // Discover merged FASTQ pairs for checking
    merged_pairs = discoverMergedFastqPairs("${params.outdir}/results/merged")
    
    log.info "Found ${merged_pairs.size()} merged FASTQ pairs for verification"
    
    // Create channel for pair checking
    check_ch = Channel.fromList(
	merged_pairs.collect { sample_id, files ->
            [sample_id, files[0], files[1]]
	}
    )
    
    // Check FASTQ pairs synchronization
    CHECK_FASTQ_PAIRS(check_ch)
    
    // Aggregate check results
    AGGREGATE_FASTQ_CHECKS(CHECK_FASTQ_PAIRS.out.check_results.map { it[1] }.collect())
    
    log.info "FASTQ pair verification completed"
    
}
