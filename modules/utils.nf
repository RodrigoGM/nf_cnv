// Discover FASTQ files in seqdata/ input directory
def discoverFastqFiles(seqdata_path, sample_pattern) {
    def samples = [:]
    
    // Find all R1 files efficiently
    def r1_files = file("${seqdata_path}/**/*_R1_*.fastq.gz").sort { it.name }
    
    r1_files.each { r1_file ->
        // Extract sample ID before Illumina conventions
        def filename = r1_file.name
        def sample_id = filename.replaceAll(/_S\d+_L\d+_R[12]_\d+\.fastq\.gz$/, '')
        
        // Apply sample pattern filter
        if (sample_id ==~ sample_pattern) {
            // Find corresponding R2 file
            def r2_file = file(r1_file.toString().replace('_R1_', '_R2_'))
            
            if (r2_file.exists()) {
                if (!samples.containsKey(sample_id)) {
                    samples[sample_id] = [r1: [], r2: []]
                }
                samples[sample_id].r1 << r1_file
                samples[sample_id].r2 << r2_file
            }
        }
    }
    
    return samples
}

// generate sample manifest
def createManifestData(samples) {
    def manifest_lines = ['sample_id,read_type,source_files,output_file']
    
    samples.each { sample_id, files ->
        def r1_sources = files.r1.collect { it.name }.join(';')
        def r2_sources = files.r2.collect { it.name }.join(';')
        
        manifest_lines << "${sample_id},R1,\"${r1_sources}\",${sample_id}_R1.fastq.gz"
        manifest_lines << "${sample_id},R2,\"${r2_sources}\",${sample_id}_R2.fastq.gz"
    }
    
    return manifest_lines.join('\n')
}


// Create Merged Fastq Pairs
def discoverMergedFastqPairs(merged_path) {
    def pairs = [:]
    
    // Find all R1 merged files
    def r1_files = file("${merged_path}/*_1.merged.fastq.gz").sort { it.name }
    
    r1_files.each { r1_file ->
        def sample_id = r1_file.name.replaceAll(/_1\.merged\.fastq\.gz$/, '')
        def r2_file = file(r1_file.toString().replace('_1.merged.fastq.gz', '_2.merged.fastq.gz'))
        
        if (r2_file.exists()) {
            pairs[sample_id] = [r1_file, r2_file]
        }
    }
    
    return pairs
}

// load barcodes
def loadBarcodes(barcode_file) {
    def barcodes = []
    def barcode_path = file(barcode_file)
    
    if (!barcode_path.exists()) {
        error "Barcode file not found: ${barcode_path}"
    }
    
    barcode_path.eachLine { line ->
        def fields = line.trim().split(/\s+/)
        if (fields.size() >= 2) {
            barcodes.add([fields[0], fields[1]])  // [barcode_name, barcode_sequence]
        }
    }
    
    return barcodes
}

