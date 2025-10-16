# Changelog

All notable changes to nf-cnv will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
- CNV Results aggregation
- Other BAM-level QC
- QC Integration with MulitQC


### Added
- GitHub repository setup with comprehensive documentation
- Example configurations and test datasets

## [0.5.0] - 2025-10-15

### Added
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
