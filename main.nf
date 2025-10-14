#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Import modules
include { CAT_FASTQ } from './modules/nf-core/cat/fastq/main'
include { discoverFastqFiles; createManifestData; loadBarcodes } from './modules/utils'
include { CHECK_FASTQ_PAIRS; AGGREGATE_FASTQ_CHECKS; BARCODE_SPLIT;
	 CHECK_CELL_FASTQ_PAIRS; AGGREGATE_CELL_CHECKS; AGGREGATE_BARCODE_STATS } from './modules/demultiplexing'
// QC Imports
include { FASTQC } from './modules/nf-core/fastqc/main'
include { FASTQSCREEN_FASTQSCREEN } from './modules/nf-core/fastqscreen/fastqscreen/main'
// include { FASTQ_SCREEN ; AGGREGATE_FASTQ_SCREEN } from './modules/qc' // (non nf-core option)

// Read Alignment
include { BWA_MEM; COLLATE_BAM; FIXMATE_BAM; SORT_BAM; MARK_DUPLICATES;
	 REMOVE_DUPLICATES; FILTER_AND_EXTRACT_READS; VALIDATE_BAM } from './modules/alignment'

// VARBIN
include { GET_BIN_COUNTS; CNV_PROFILE; AGGREGATE_CNV_RESULTS } from './modules/varbin'
	
// BAM QC 
include { QUALIMAP_BAMQC } from './modules/nf-core/qualimap/bamqc/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from './modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTINSERTSIZEMETRICS } from './modules/nf-core/picard/collectinsertsizemetrics/main'
include { SAMTOOLS_FLAGSTAT } from './modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_STATS } from './modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from './modules/nf-core/samtools/idxstats/main'


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Global Parameters
params.seqdata = "${launchDir}/seqdata"
params.outdir = "${launchDir}"
params.sample_pattern = ".*"
params.barcode_file = "${launchDir}/nf_cnv/assets/barcode.192.txt"


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
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


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
// Main workflow
workflow {
    log.info " ███▄    █   █████▒      ▄████▄   ███▄    █ ██▒   █▓ "
    log.info " ██ ▀█   █ ▓██   ▒      ▒██▀ ▀█   ██ ▀█   █▓██░   █▒ "
    log.info "▓██  ▀█ ██▒▒████ ░      ▒▓█    ▄ ▓██  ▀█ ██▒▓██  █▒░ "
    log.info "▓██▒  ▐▌██▒░▓█▒  ░      ▒▓▓▄ ▄██▒▓██▒  ▐▌██▒ ▒██ █░░ "
    log.info "▒██░   ▓██░░▒█░         ▒ ▓███▀ ░▒██░   ▓██░  ▒▀█░   "
    log.info "░ ▒░   ▒ ▒  ▒ ░         ░ ░▒ ▒  ░░ ▒░   ▒ ▒   ░ ▐░   "
    log.info "░ ░░   ░ ▒░ ░             ░  ▒   ░ ░░   ░ ▒░  ░ ░░   "
    log.info "   ░   ░ ░  ░ ░         ░           ░   ░ ░     ░░   "
    log.info "         ░              ░ ░               ░      ░   "
    log.info "                        ░                      ░     "
    log.info "."
    log.info "."
    log.info "SC-CNV Processing Pipeline"
    log.info "Using genome: ${params.genome}"
    
    // Validate genome parameter
    def valid_genomes = ['hsa37', 'hsa38', 'pdx37', 'mm10', 'mm38', 'mm39']
    if (!(params.genome in valid_genomes)) {
        error "Invalid genome: ${params.genome}. Valid options: ${valid_genomes.join(', ')}"
    }
    
    // Verify genome files exist
    def bwa_index = params.genome_indices[params.genome]
    def ref_fasta = params.genome_fasta[params.genome]
    
    if (!file(bwa_index).exists()) {
        error "BWA index not found: ${bwa_index}"
    }
    if (!file(ref_fasta).exists()) {
        error "Reference FASTA not found: ${ref_fasta}"
    }
    
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
    
    log.info "Demultiplexing complete: ${samples.size()} samples, ${barcodes.size()} barcodes per sample"

    // STEP 5: CELL-LEVEL QUALITY CONTROL
    // Run FastQC on demultiplexed cell FASTQ files
    cell_fastq_for_qc = BARCODE_SPLIT.out.cell_fastqs.map { cell_id, r1, r2 -> 
	[[id: cell_id, single_end: false], [r1, r2]]
    }
    FASTQC(cell_fastq_for_qc)
    
    log.info "FastQC completed on demultiplexed cells"

    // FastQ Screen QC (on demultiplexted FASTQs)
    if (params.run_fastq_screen) {
	// Prepare input channel from split cell fastqs
	cell_fastqs_for_screen = BARCODE_SPLIT.out.cell_fastqs
            .map { cell_id, r1, r2 -> [[id: cell_id, single_end: false], [r1, r2]] }
	
	// Create empty database channel if no config provided
	fastq_screen_db = params.fastq_screen_config ? 
            Channel.fromPath(params.fastq_screen_config) : 
            Channel.empty()
	
	FASTQSCREEN_FASTQSCREEN(
            cell_fastqs_for_screen,
            fastq_screen_db
	)
    }
    
    // STEP 6: SEQUENCE ALIGNMENT
    // Run BWA-MEM alignment on cell FASTQ files
    alignment_input = BARCODE_SPLIT.out.cell_fastqs
    BWA_MEM(alignment_input)

    // Validate BAM files before processing
    VALIDATE_BAM(BWA_MEM.out.unsorted_bam)

    // Process BAM files through standard pipeline
    COLLATE_BAM(VALIDATE_BAM.out.validated_bam)
    FIXMATE_BAM(COLLATE_BAM.out.collated_bam)
    SORT_BAM(FIXMATE_BAM.out.fixmate_bam)

    // Parallel processing for marked and deduplicated BAMs
    MARK_DUPLICATES(SORT_BAM.out.sorted_bam)
    REMOVE_DUPLICATES(MARK_DUPLICATES.out.marked_bam)

    log.info "Alignment workflow completed for genome: ${params.genome}"
    
    // STEP 7: POST-PROCESS DEDUPLICATED BAMS
    // Filter for unique reads and extract forward/reverse strand
    strand_suffix = params.use_reverse_reads ? 'RV' : 'FW'
    FILTER_AND_EXTRACT_READS(REMOVE_DUPLICATES.out.dedup_bam)
    
    log.info "Post-alignment filtering complete: extracted ${strand_suffix} reads for single-end analysis"

    // BAM Quality Control - CORRECTED VERSION
    if (params.run_bam_qc) {
	// Use deduplicated BAMs for QC
	qc_bams = REMOVE_DUPLICATES.out.dedup_bam
	
	// QualiMap BAM QC
	QUALIMAP_BAMQC(
            qc_bams.map { cell_id, bam, bai -> [cell_id, bam] },
            []  // No GFF file
	)
	
	// Picard metrics
	PICARD_COLLECTMULTIPLEMETRICS(
            qc_bams.map { cell_id, bam, bai -> [cell_id, bam] },
            [], // No reference fasta needed for basic metrics
            []  // No reference dict needed
	)
	
	PICARD_COLLECTINSERTSIZEMETRICS(
            qc_bams.map { cell_id, bam, bai -> [cell_id, bam] }
	)
	
	// Samtools metrics
	SAMTOOLS_FLAGSTAT(
            qc_bams.map { cell_id, bam, bai -> [cell_id, bam, bai] }
	)
	
	SAMTOOLS_STATS(
            qc_bams.map { cell_id, bam, bai -> [cell_id, bam] },
            [] // No reference fasta
	)
	
	SAMTOOLS_IDXSTATS(
            qc_bams.map { cell_id, bam, bai -> [cell_id, bam, bai] }
	)
    }
    
    // Varbin CNV workflow
    if (params.run_cnv_analysis &&
	params.genome in params.varbin_supported_genomes) {

	// Varbin CNV
	// Create resolution channel
	resolution_ch = Channel.from(params.cnv_resolutions)
	
	// Combine strand BAMs with resolutions
	cnv_input_ch = FILTER_AND_EXTRACT_READS.out.strand_bam
            .combine(resolution_ch)
	
	// Run bin counting
	GET_BIN_COUNTS(cnv_input_ch)
	
	// Run CNV profiling
	CNV_PROFILE(GET_BIN_COUNTS.out.bin_counts)
    
	// Aggregate results across all cells and resolutions
	all_ploidy_results = CNV_PROFILE.out.ploidy_results
            .map { cell_id, resolution, ploidy_file -> ploidy_file }
            .collect()
    
	AGGREGATE_CNV_RESULTS(all_ploidy_results)
    }
    
}


//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
def helpMessage() {
    log.info"""
    =======================================
    SC-CNV Nextflow Pipeline
    =======================================

    Usage:
    nextflow run nf_cnv/main.nf -c config/nextflow.config [options]

    Options:
      --seqdata         Path to sequencing data directory (default: ${params.seqdata})
      --outdir          Output directory (default: ${params.outdir})
      --sample_pattern  Sample filtering regex (default: ${params.sample_pattern})
      --genome          Reference genome [hsa37|hsa38|pdx37|mm10|mm38|mm39] (default: ${params.genome})
	... [nextflow args]
      --help            Display this help message

    Examples:
      # Human genome GRCh38
      nextflow run nf_cnv/main.nf -c config/nextflow.config --genome hsa38
      
      # Mouse genome mm39
      nextflow run nf_cnv/main.nf -c config/nextflow.config --genome mm39
      
      # Filter samples and use specific genome
      nextflow run nf_cnv/main.nf -c config/nextflow.config --genome hsa37 --sample_pattern "LS0775.*"
    """.stripIndent()
}
