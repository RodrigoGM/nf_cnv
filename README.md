# nf-cnv: Single Cell Copy Number Variation Analysis Pipeline

[![nf-cnv v0.7.3](https://img.shields.io/badge/version-0.7.3-D40000.svg)](https://github.com/RodrigoGM/nf_cnv/)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![License: BSD-3](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/bsd-3-clause)


A comprehensive Nextflow pipeline for single-cell copy number variation (CNV) analysis from multiplexed FASTQ files.

## Pipeline Overview

nf-cnv processes multiplexed single-cell sequencing data through the complete CNV analysis workflow:

FASTQ Files → Demultiplexing → Quality Control → Alignment → CNV Calling → Summary Statistics


### Key Features

- **Automated demultiplexing** of multiplexed single-cell FASTQ files using custom barcodes
- **Comprehensive quality control** with FastQC, FastQ Screen, and BAM-level QC
- **Robust alignment pipeline** using BWA-MEM with duplicate removal and filtering
- **Multi-resolution CNV analysis** using varbin approach (5k, 20k, 50k bins)
- **Strand-specific analysis** with forward/reverse read extraction
- **Automated statistics** and contamination reporting
- **Scalable execution** on SLURM clusters with resource optimization

## Quick Start

### 1. Installation

```bash
# Clone the repository
git clone https://github.com/RodrigoGM/nf_cnv.git
cd nf_cnv

# Set up conda environment
conda env create -f environment.yml
conda activate single-cell-cnv

```

### 2. Initiate a project
```bash
# Create project directory
mkdir my_project && cd my_project

# Initialize project structure
../nf_cnv/bin/init.sh
```

### 3. Configure Your Analysis
Edit the generated configuration files:
* `config/nextflow.config` - Pipeline parameters
* `config/samples.csv` - Sample metada

### 4. Run the pipeline


``` bash
# Test run
./run_pipeline.sh --dry-run

# Full execution
./run_pipeline.sh
# or submit to SLURM
sbatch run_pipeline.sh
```


## Pipeline Workflow

### 1. FASTQ Processing
 **Concatenation** : Merge FASTQ files from multiple sequencing runs
 **Validation** : Verify paired-end synchronization
 **Demultiplexing** : Split by single-cell barcodes (192 barcodes default)

### 2. Quality Control
 **FastQC** : Per-cell quality metrics
 **FastQ Screen** : Contamination detection
 **BAM QC** : Alignment statistics with Picard and Samtools

### 3. Sequence Alignment
 **BWA-MEM** : Paired-end alignment to reference genome
 **Processing** : Collation, fixmate, sorting, duplicate marking/removal
 **Filtering** : Strand-specific read extraction with chromosome filtering

### 4. CNV Analysis
 **Bin Counting** : Read depth calculation across genomic bins
 **CNV Processing** : Using variable length binning (Varbin)
 **Ploidy Estimation** : Ploidy estimation using leaset squares regression
 **Multi-resolution** : Analysis at 5kb, 20kb, and 50kb resolution
 **Statistics** : Comprehensive summary and quality metrics

## Directory Structure

``` bash
project/
├── config/
│   ├── nextflow.config     # Main configuration
│   ├── samples.csv         # Sample metadata
│   └── fastq_screen.conf   # FastQ Screen databases
├── seqdata/                # Raw FASTQ files
├── results/                # Pipeline outputs
│   ├── wsplit/            # Demultiplexed cells
│   ├── bwa_out_md/        # Marked Duplicates BAMs
│   ├── bwa_out_dd/        # DeDuplicated BAMs
│   ├── bwa_out_fw/        # Read 1 BAMs Q30, Primary Genome
|	├── varbin5k/          # CNV results (5kb)
|	│   ├── counts/             # *.bin.counts.bed
|	│   ├── stats/              # *.bin.counts.stats.bed
|	│   ├── seg_long/           # *_seg.txt
|	│   ├── seg_short/          # *_short_seg.txt
|	│   └── ploidy/             # *.quantal.ploidy.txt
|	├── varbin20k/         # CNV results (20kb)
|	│   └── ...
|	└── varbin50k/         # CNV results (50kb)
|	|   └── ...
│   ├── cnv_summary/       # Results and Statistics Summaries
│   └── qc/                # Quality control reports
├── logs/                  # Execution logs and reports
└── nf_cnv/                # Pipeline code (symlink)
```

## Configuration
| Parameter         | Description                                           | Default |
|:------------------|:------------------------------------------------------|---------|
| genome            | Reference genome [hsa37,hsa38,mm10,mm39, pdx37,pdx38] | hsa38   |
| use_reverse_reads | Extract reverse reads (R2) instead of forward         | FALSE   |
| min_ploidy        | Minimum ploidy for CNV calling                        | 1.8     |
| max_ploidy        | Maximum ploidy for CNV calling                        | 4.8     |
| run_cnv_analysis  | Enable CNV analysis                                   | TRUE    |
| run_bam_qc        | Enable BAM quality control                            | TRUE    |


## Genome Support
 **Human** : hsa37 (GRCh37), hsa38 (GRCh38)
 *Mouse* : mm10, mm39 (pending cnv)
 *PDX models* : pdx37, pdx38 (human-mouse hybrid genome) (pending cnv)
 

## Usage Examples
### Basic Human Analysis

``` bash
./run_pipeline.sh --genome hsa38

```
### Mouse Analysis

``` bash
./run_pipeline.sh --genome mm39
```

### Reverse Read Analysis

``` bash
./run_pipeline.sh --use-reverse-reads
```

### Custom Ploidy Range
./run_pipeline.sh --min-ploidy 1.5 --max-ploidy 5.0

### Specific Sample Pattern
./run_pipeline.sh --sample-pattern "LS0775.*"

## Output Files

### CNV Results
 *.quantal.ploidy.txt - Ploidy estimates per cell
 *.bin.counts.bed - Read counts per genomic bin
 *_seg.txt - Segmentation results
 *_short_seg.txt - Condensed segmentation

### Quality Control
 fastq_screen_summary.txt - Contamination analysis
 *_flagstat - Alignment statistics
 *.insert_size_metrics - Insert size distributions
### Summary Statistics
 ploidy_results_combined.txt - All ploidy results
 ploidy_summary_all.txt - Statistical summaries
 contamination_report.txt - QC flagged cells

## Requirements

### Software Dependencies
 Nextflow ≥23.04.0
 Conda/Mamba
 SLURM (for cluster execution)

### Reference Data
 BWA indices for supported genomes
 cna_utils reference files
 FastQ Screen databases (optional)

### Computational Resources
 Memory : 8-64GB per process
 CPUs : 1-8 cores per process
 Storage : 1.3TB per 1000 cells (@ 2.4M reads / cell)
 Time : 4-24 hours for typical datasets

## Troubleshooting
### Common Issues

#### Pipeline fails to start
``` bash
# Check configuration
nextflow run nf_cnv/main.nf -c config/nextflow.config --dry-run

# Verify sample sheet
head config/samples.csv
```

#### SLURM job failures
* Check account and partition settings
* Verify resource allocations in conf/modules.config
* Review SLURM logs in logs/

#### Missing output files
* Check process work directories: work/*/
* Verify publishDir patterns in configuration
* Use --resume to restart from last successful step
#### Getting Help
1. Check execution reports in logs/reports/
2. Review .nextflow.log for detailed errors
3. Open issues on GitHub

## Citation

``` bash
nf-cnv: A Nextflow pipeline for single-cell copy number variation analysis
https://github.com/RodrigoGM/nf_cnv
```
## License
This project is licensed under the BSD-3 Clause License - see the LICENSE file for details.

## Contributing
Submit a documeted pull request

## Acknowledgments
* _Nextflow_ for workflow management
* _nf-core_ for module templates
* [cna_utils](https://github.com/rishvanth-kp/cna_utils) for CNV analysis tools

