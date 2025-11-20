# Changelog

All notable changes to nf-cnv will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- **QDNAseq Analysis**: 
  - include 100 and 150 mappable Kb variable length bins for mm39 and hsa38


- **Preseq Analysis**: Library complexity estimation suite
  - C-curve complexity analysis
  - LC-extrap library complexity extrapolation
  - GC-extrap genomic coverage extrapolation
  - Separate output directories for each metric

- Input verification for `GET_BIN_COUNTS`

- Conditional VARBIN outputs

## [0.7.4] - 2025-11-20

### Changed
- Conditional FASTQ aligment, only cells with greater than params.minimum_read_count are aligned

### Fixed
- Phenotype template file (`cell_phenotype.txt`) now includes *only* cells actually passing read-count thresholds and processed/aligned, ensuring downstream consistency and avoiding empty rows.
- Robust join/filter for sequence alignment using barcode stats, protecting pipeline against sample filtering mismatches.
- main.nf code formatting


## [0.7.3] - 2025-10-28

### Changed
- **Conditional Outputs**: QC, and intermediary BAM files are no longer published by default.  Varbin process outputs are now routed to new, organized subdirectories, however plans are being made to conditionally publish.

- **Slurm partition usage**: Most processes are now joinly submitted two paritions (a standard and a preemptable by default) for improved resource cost-efficiency without affecting core results.

- **Process cleanup and staging**: Enhanced intermediate file cleanup and staging delays for better stability and performance on large datasets and NFS-backed clusters.

### Fixed
- Small bugs and redundancies in output and resource management.


## [0.7.2] - 2025-10-28

### Fixed
- **NFS File Staging Reliability**: Enhanced robustness for large-scale datasets
  - Increased file staging timeout to 30 minutes
  - Added 5-retry mechanism with exponential backoff for file transfers
  - Implemented sync delays (3-5s) before process execution in I/O-heavy steps
  - Added input validation in GET_BIN_COUNTS to detect missing BAM files
  
- **Error Handling Strategy**: Unified retry logic
  - Ignore QC non-recoverable errors to maintain pipeline flow
  - Prevent premature pipeline termination on isolated failures

- **File Porter Configuration**: Optimized for NFS environments
  - maxRetries increased from 3 to 5
  - pollInterval adjusted to 15 seconds to reduce NFS load
  - maxTransfers set to 50 to prevent filesystem overwhelm

### Changed
- Simplified BARCODE_SPLIT dummy read generation logic
- Increase staging timing  `CAT_FASTQ`, `FILTER_AND_EXTRACT_READS`, and `GET_BIN_COUNTS`
- Consolidated error handling across all processes
- Enhanced beforeScript, and afterScript delays for large file operations


## [0.7.1] - 2025-10-22

### Added
- **MultiQC Integration** : 
  - create process dependencies so it runs after all QC processes are complete

## [0.7.0] - 2025-10-18

### Added
- **MultiQC Integration**: Comprehensive QC report aggregation with nf-core module
  - Custom report naming with project-specific titles
  - Output to cnv_summary directory with uppercase naming convention

- **Cell Phenotype Templates**: Automated metadata template generation
  - Pre-filled with cellID, ploidy, and error from 20k resolution
  - Template columns for experimental annotation
  - Published to cnv_summary directory

- **Structured Varbin Outputs**: Organized CNV result directories
  - counts/ and stats/ subdirectories for bin counts
  - ploidy/, seg_long/, seg_short/ subdirectories for CNV calls
  - Improved result organization and accessibility

### Enhanced
- **Configuration Consolidation**: Unified resource management
  - Consolidated process resources into modules.config
  - Clearer separation: params in base.config, resources in modules.config

- **Resource Optimization**: Improved memory and CPU allocation
  - MultiQC: 32GB for large dataset processing
  - AGGREGATE_SEG_MATRICES: Configurable memory (default 32GB)
  - FastQ Screen: 64GB for large genome indices

- **Work Directory Cleanup**: Automatic intermediate file removal
  - Reduced work directory footprint
  - Cleanup of merged FASTQs, intermediate BAM, and temporary files
  - Configurable cleanup policies per process

### Fixed
- MultiQC output pattern matching for custom filenames
- BAM QC module input channel formatting for nf-core compatibility
- FastQ Screen error handling for empty files
- Memory allocation overrides in SLURM configuration
- Process resource hierarchy and precedence rules

### Performance
- Optimized matrix aggregation with memory-efficient bash implementation
- Automatic cleanup reduces storage requirements by ~70%
- Improved SLURM resource utilization with proper queue management

### Infrastructure
- Enhanced init.sh with improved project setup workflow
- Updated configuration hierarchy documentation
- Clarified params vs process config separation


## [0.6.0] - 2025-10-17

### Added
- **Result Organization**: Structured output directories for varbin results
  - `counts/` and `stats/` subdirectories for bin count files
  - `ploidy/`, `seg_long/`, `seg_short/` subdirectories for CNV outputs

- **Matrix Aggregation**: Memory-efficient seg.txt matrix generation
  - Bash-based implementation for handling 10K+ cell datasets
  - Separate matrices for counts, ratios, and segmentation results

- **Cell Phenotype Templates**: Automated phenotype table generation
  - Pre-filled with cellID, ploidy, and error from 20k resolution
  - Template columns for experimental metadata annotation

- **Enhanced QC Workflows**: Improved quality control processes
  - FastQ Screen integration with empty file handling
  - Comprehensive BAM QC with Picard and Samtools modules
  - Proper nf-core module integration and error strategies

### Enhanced
- **Resource Management**: Optimized memory allocation for large datasets
- **Error Handling**: Robust handling of empty FASTQ and BAM files
- **Process Organization**: Structured publishDir patterns for better file management
- **Chromosome Filtering**: Dynamic chromosome selection for strand-specific analysis

### Fixed
- FastQ Screen failures on zero-read files using errorStrategy ignore
- BAM QC module input channel formatting for nf-core compatibility
- Memory exhaustion in seg matrix aggregation for large cell numbers
- Cell phenotype template generation with proper 20k resolution filtering

### Performance
- Streamlined matrix operations for datasets with 1000+ cells
- Reduced memory footprint for aggregation processes
- Improved SLURM resource utilization and job scheduling

### Infrastructure
- Updated module configurations with appropriate resource allocations
- Enhanced debugging capabilities with process execution logging
- Comprehensive documentation updates with usage examples


## [0.5.0] - 2025-10-15

### Added
- GitHub repository setup with comprehensive documentation
- Example configurations and test datasets
- Complete single-cell CNV processing pipeline
- Automated barcode demultiplexing with 192-barcode support
- Multi-resolution CNV calling (5k, 20k, 50k bins)
- Comprehensive quality control suite
- Custom FastQ Screen implementation
- BAM-level QC with Picard and Samtools modules
- Strand-specific read extraction and filtering
- Automated statistics and summary reporting
- SLURM cluster optimization with resource management
- Interactive project initialization script
- Python-based ploidy statistics analysis
- Support for multiple reference genomes (human, mouse, PDX)

### Pipeline Components
- **FASTQ Processing**: Concatenation, validation, demultiplexing
- **Quality Control**: FastQC, FastQ Screen, contamination detection
- **Alignment**: BWA-MEM with full BAM processing pipeline
- **CNV Analysis**: varbin-based copy number calling
- **Reporting**: HTML reports, statistics, and visualization

### Supported Genomes
- Human: hsa37 (GRCh37), hsa38 (GRCh38)
- Mouse: mm10, mm39
- PDX: pdx37 (hg37+mm10), pdx38(hg38+mm39)

### Features
- Configurable ploidy ranges (default: 1.8-4.8)
- Forward/reverse read selection for single-end analysis
- Chromosome filtering for targeted analysis
- Comprehensive error handling and resume capability
- Scalable execution with resource optimization

## [0.4.0] - 2025-10-10

### Added
- Core alignment pipeline with BWA-MEM
- Basic CNV calling functionality
- SLURM executor configuration
- Initial barcode splitting implementation

### Fixed
- Channel mapping issues in workflow
- Process input/output tuple mismatches
- FASTQ file discovery and validation

## [0.3.0] - 2025-10-05

### Added
- Initial Nextflow DSL2 pipeline structure
- nf-core module integration
- Basic FASTQ processing workflows

### Infrastructure
- Conda environment setup
- Base configuration templates
- Module organization structure

## [0.1.0] - 2025-09-01

### Added
- Project
