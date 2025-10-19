# Changelog

All notable changes to nf-cnv will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
- QC Integration with MulitQC


### Added

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
