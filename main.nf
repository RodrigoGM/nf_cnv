#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { CAT_FASTQ } from './modules/nf-core/cat/fastq/main'
include { discoverFastqFiles; createManifestData; loadBarcodes } from './modules/utils'
include { CHECK_FASTQ_PAIRS; AGGREGATE_FASTQ_CHECKS; BARCODE_SPLIT; CHECK_CELL_FASTQ_PAIRS; AGGREGATE_CELL_CHECKS; AGGREGATE_BARCODE_STATS } from './modules/demultiplexing'

// Global Parameters
params.seqdata = "${launchDir}/seqdata"
params.outdir = "${launchDir}"
params.sample_pattern = ".*"
params.barcode_file = "${launchDir}/nf_cnv/assets/barcode.192.txt"

// Process for creating FASTQ manifest
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

// Main workflow
workflow {
    log.info "SC-CNV Processing Pipeline"
    
    // STEP 1: FASTQ CONCATENATION
    // Discover and group FASTQ files by sample ID
    samples = discoverFastqFiles(params.seqdata, params.sample_pattern)
    log.info "Found ${samples.size()} samples matching pattern '${params.sample_pattern}'"
    
    if (samples.isEmpty()) {
        error "No samples found matching pattern '${params.sample_pattern}' in ${params.seqdata}"
    }
    
    // Create channel for nf-core CAT_FASTQ module
    sample_ch = Channel.fromList(
        samples.collect { sample_id, files ->
            def meta = [id: sample_id, single_end: false]
            def reads = [files.r1, files.r2].flatten().sort { it.name }
            [meta, reads]
        }
    )
    
    // Concatenate FASTQ files using nf-core module
    CAT_FASTQ(sample_ch)
    
    // Create concatenation manifest
    manifest_data = createManifestData(samples)
    CREATE_FASTQ_MANIFEST(Channel.value(manifest_data))
    
    // STEP 2: FASTQ PAIR VERIFICATION
    // Verify merged FASTQ pairs are synchronized
    merged_pairs_ch = CAT_FASTQ.out.reads.map { meta, reads ->
        [meta.id, reads[0], reads[1]]
    }
    
    CHECK_FASTQ_PAIRS(merged_pairs_ch)
    AGGREGATE_FASTQ_CHECKS(CHECK_FASTQ_PAIRS.out.check_results.map { it[1] }.collect())
    log.info "FASTQ pair verification completed"
    
    // STEP 3: BARCODE DEMULTIPLEXING
    // Load barcodes for single-cell demultiplexing
    barcodes = loadBarcodes(params.barcode_file)
    log.info "Loaded ${barcodes.size()} barcodes for demultiplexing"
    
    // Create sample-barcode combinations for parallel processing
    barcode_ch = CAT_FASTQ.out.reads
        .map { meta, reads -> [meta.id, reads[0], reads[1]] }
        .combine(Channel.fromList(barcodes))
        .map { sample_id, r1, r2, barcode_name, barcode_seq ->
            [sample_id, r1, r2, barcode_name, barcode_seq]
        }
    
    // Execute barcode splitting in parallel
    BARCODE_SPLIT(barcode_ch)
    
    // STEP 4: CELL-LEVEL QC
    // Verify demultiplexed cell FASTQ pairs
    CHECK_CELL_FASTQ_PAIRS(BARCODE_SPLIT.out.cell_fastqs)
    
    // Aggregate all results
    AGGREGATE_CELL_CHECKS(CHECK_CELL_FASTQ_PAIRS.out.collect())
    AGGREGATE_BARCODE_STATS(BARCODE_SPLIT.out.stats.map { it[2] }.collect())
    
    log.info "Pipeline completed: ${samples.size()} samples, ${barcodes.size()} barcodes per sample"
}
