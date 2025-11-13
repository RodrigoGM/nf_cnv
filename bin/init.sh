#!/bin/bash
# nf_cnv/init.sh
# Interactive script to initialize a project for the SC-CNV Nextflow pipeline
# Generates config files, run script, and project documentation

# Set text colors
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Put colors into an array (exclude NC)
COLORS=("$RED" "$GREEN" "$BLUE" "$YELLOW" "$PURPLE" "$CYAN")

# Pick a random color from the array
RANDOM_COLOR=${COLORS[$RANDOM % ${#COLORS[@]}]}

# ASCII art banner
echo -e "${RANDOM_COLOR}"
cat << "EOF"
 ███▄    █   █████▒      ▄████▄   ███▄    █ ██▒   █▓
 ██ ▀█   █ ▓██   ▒      ▒██▀ ▀█   ██ ▀█   █▓██░   █▒
▓██  ▀█ ██▒▒████ ░      ▒▓█    ▄ ▓██  ▀█ ██▒▓██  █▒░
▓██▒  ▐▌██▒░▓█▒  ░      ▒▓▓▄ ▄██▒▓██▒  ▐▌██▒ ▒██ █░░
▒██░   ▓██░░▒█░         ▒ ▓███▀ ░▒██░   ▓██░  ▒▀█░  
░ ▒░   ▒ ▒  ▒ ░         ░ ░▒ ▒  ░░ ▒░   ▒ ▒   ░ ▐░  
░ ░░   ░ ▒░ ░             ░  ▒   ░ ░░   ░ ▒░  ░ ░░  
   ░   ░ ░  ░ ░         ░           ░   ░ ░     ░░  
         ░              ░ ░               ░      ░  
                        ░                      ░   
EOF
echo -e "${NC}"
echo -e "${GREEN}============================================${NC}"
echo -e "${GREEN}Single Cell CNV Pipeline Project Setup${NC}"
echo -e "${GREEN}============================================${NC}"

# Get the directory where init.sh is located (should be nf_cnv/)
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PIPELINE_DIR=$(dirname $SCRIPT_DIR)

# Determine project directory (parent of nf_cnv)
PROJECT_DIR=$(dirname $(dirname $SCRIPT_DIR))
echo -e "${BLUE}Pipeline directory:${NC} $PIPELINE_DIR"
echo -e "${BLUE}Project directory:${NC} $PROJECT_DIR"

# Extract project name from directory
PROJECT_NAME=$(basename $PROJECT_DIR)
echo -e "\n${PURPLE}Project Information:${NC}"
read -p "Enter project name (default: $PROJECT_NAME): " PROJECT_NAME_INPUT
if [ -n "$PROJECT_NAME_INPUT" ]; then
    PROJECT_NAME="$PROJECT_NAME_INPUT"
fi

read -p "Enter project description: " PROJECT_DESCRIPTION
read -p "Enter your name/organization: " PROJECT_AUTHOR
read -p "Enter contact email: " PROJECT_EMAIL

# Verify seqdata directory exists
SEQDATA_DIR="$PROJECT_DIR/seqdata"
if [ ! -d "$SEQDATA_DIR" ]; then
    echo -e "${YELLOW}Warning: seqdata directory not found at $SEQDATA_DIR${NC}"
    read -p "Enter path to seqdata directory (or press Enter to create it): " SEQDATA_INPUT
    if [ -n "$SEQDATA_INPUT" ]; then
        SEQDATA_DIR="$SEQDATA_INPUT"
        if [ ! -d "$SEQDATA_DIR" ]; then
            mkdir -p "$SEQDATA_DIR"
            echo -e "${GREEN}Created directory: $SEQDATA_DIR${NC}"
        fi
    else
        mkdir -p "$SEQDATA_DIR"
        echo -e "${GREEN}Created directory: $SEQDATA_DIR${NC}"
    fi
fi

# Create necessary directories
CONFIG_DIR="$PROJECT_DIR/config"

echo -e "\n${BLUE}Creating project structure...${NC}"
mkdir -p "$CONFIG_DIR"

# Function to extract sample IDs from FASTQ files
function find_sample_ids() {
    echo -e "\n${BLUE}Searching for FASTQ files in $SEQDATA_DIR...${NC}"
    
    # Look for fastq files recursively
    FASTQ_FILES=$(find "$SEQDATA_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | sort)
    
    if [ -z "$FASTQ_FILES" ]; then
        echo -e "${YELLOW}No FASTQ files found. Sample sheet will be created with example entries.${NC}"
        return 1
    fi
    
    echo -e "${GREEN}Found $(echo "$FASTQ_FILES" | wc -l) FASTQ files.${NC}"
    
    # Extract sample IDs from filenames using multiple patterns
    SAMPLE_IDS=""
    
    # Pattern 1: SampleName_S1_L001_R1_001.fastq.gz
    PATTERN1=$(echo "$FASTQ_FILES" | grep -o -E "[A-Za-z0-9_-]+_S[0-9]+_L[0-9]+_R[12]_[0-9]+" | sed 's/_S[0-9]*_L[0-9]*_R[12]_[0-9]*//' | sort | uniq)
    
    # Pattern 2: SampleName_R1.fastq.gz
    PATTERN2=$(echo "$FASTQ_FILES" | grep -o -E "[A-Za-z0-9_-]+_R[12]" | sed 's/_R[12]//' | sort | uniq)
    
    # Pattern 3: SampleName.fastq.gz
    PATTERN3=$(echo "$FASTQ_FILES" | xargs -n 1 basename | sed 's/\.\(fastq\|fq\)\.gz$//' | grep -v "_R[12]" | sort | uniq)
    
    # Choose the best pattern
    if [ -n "$PATTERN1" ]; then
        SAMPLE_IDS="$PATTERN1"
        echo -e "${GREEN}Using Illumina naming pattern${NC}"
    elif [ -n "$PATTERN2" ]; then
        SAMPLE_IDS="$PATTERN2"
        echo -e "${GREEN}Using simple R1/R2 pattern${NC}"
    elif [ -n "$PATTERN3" ]; then
        SAMPLE_IDS="$PATTERN3"
        echo -e "${GREEN}Using basic filename pattern${NC}"
    else
        echo -e "${RED}Failed to automatically extract sample IDs.${NC}"
        return 1
    fi
    
    echo -e "${GREEN}Found $(echo "$SAMPLE_IDS" | wc -l) unique sample IDs:${NC}"
    echo "$SAMPLE_IDS" | sed 's/^/  - /'
    return 0
}

# Gather SLURM configuration details
echo -e "\n${PURPLE}SLURM Configuration:${NC}"
read -p "SLURM account name (default: singers): " SLURM_ACCOUNT
SLURM_ACCOUNT=${SLURM_ACCOUNT:-singers}

read -p "SLURM queue/partition (default: cpu): " SLURM_QUEUE
SLURM_QUEUE=${SLURM_QUEUE:-cpu}

read -p "Maximum concurrent jobs (default: 50): " MAX_JOBS
MAX_JOBS=${MAX_JOBS:-50}

read -p "Submit rate limit (default: '10 sec'): " SUBMIT_RATE
SUBMIT_RATE=${SUBMIT_RATE:-'10 sec'}

# CPU/Memory configuration
echo -e "\n${PURPLE}Resource Configuration:${NC}"
read -p "CPUs for concat process (default: 2): " CONCAT_CPUS
CONCAT_CPUS=${CONCAT_CPUS:-2}

read -p "Memory for concat process in GB (default: 16): " CONCAT_MEM
CONCAT_MEM=${CONCAT_MEM:-16}

read -p "CPUs for barcode splitting (default: 4): " SPLIT_CPUS
SPLIT_CPUS=${SPLIT_CPUS:-4}

read -p "Memory for barcode splitting in GB (default: 8): " SPLIT_MEM
SPLIT_MEM=${SPLIT_MEM:-8}

# Barcode Configuration
echo -e "\n${PURPLE}Barcode Configuration:${NC}"
read -p "Barcode file path (default: \${launchDir}/nf_cnv/assets/barcode.192.txt): " BARCODE_FILE
BARCODE_FILE=${BARCODE_FILE:-"\${launchDir}/nf_cnv/assets/barcode.192.txt"}

read -p "Maximum barcode mismatches allowed (default: 1): " BARCODE_MISMATCHES
BARCODE_MISMATCHES=${BARCODE_MISMATCHES:-1}

read -p "Barcode search length (default: 80): " BARCODE_SEARCH_LEN
BARCODE_SEARCH_LEN=${BARCODE_SEARCH_LEN:-80}

read -p "Enable reverse complement matching? (y/n, default: y): " BARCODE_RC_INPUT
BARCODE_RC="true"
if [[ "$BARCODE_RC_INPUT" == "n" ]]; then
    BARCODE_RC="false"
fi

read -p "Trim barcode from reads? (y/n, default: n): " BARCODE_TRIM_INPUT
BARCODE_TRIM="false"
if [[ "$BARCODE_TRIM_INPUT" == "y" ]]; then
    BARCODE_TRIM="true"
fi

echo -e "\n${BLUE}Generating configuration files...${NC}"

# Create nextflow.config
cat > "$CONFIG_DIR/nextflow.config" << EOL
/*
 * Nextflow configuration for $PROJECT_NAME
 * Generated on: $(date)
 * Author: $PROJECT_AUTHOR
 */

// Include base configurations
includeConfig "\${launchDir}/nf_cnv/conf/base.config"
includeConfig "\${launchDir}/nf_cnv/conf/modules.config"
includeConfig "\${launchDir}/nf_cnv/conf/slurm.config"

// Project-specific configuration
// Set work directory for this project
// workDir = "${launchDir}/work"


params {
    // Project Information
    project = "$PROJECT_NAME"
    description = "$PROJECT_DESCRIPTION"
    genome = "hsa38" // select genome [hsa37, mm10, mm39, pdx37]
    primary_chromosomes = "{1..22} X Y"

    // seg aggregation override
    aggregate_seg_memory = '64.GB'    
    
    // Output Options
    publish_dir_mode = 'copy'
    save_intermediate = false
}

// Process Configuration
process {
    executor = 'slurm'
    queue = { task.attempt == 1 ? 'preemptable' :'$SLURM_QUEUE' }
    clusterOptions = '--account=$SLURM_ACCOUNT'
    
    // Error handling
    errorStrategy = {
        task.exitStatus in [143, 137, 271] ? 'retry' : 'ignore'
        task.attempt <= 10 ? 'retry' : 'finish' }
    maxRetries = 10

    // Critical/fast processes always on cpu
    withName: 'AGGREGATE_*|CREATE_*|MULTIQC' {
        queue = '$SLURM_QUEUE'
        maxRetries = 10
    }

    // Long processes: split load
    withName: 'BARCODE_SPLIT|BWA_MEM|MARK_DUPLICATES|REMOVE_DUPLICATES|GET_BIN_COUNTS|CNV_PROFILE' {
        queue = { 
            // 50/50 split based on task hash
            Math.abs(task.hash.toString().hashCode()) % 2 == 0 ? 
                '$SLURM_QUEUE' : 'preemptable' 
        }
        clusterOptions = { task.attempt > 6 ? 
            "--partition='$SLURM_QUEUE'" : ""  // Retries always go to $SLURM_QUEUE
        }
        maxRetries = 10
    }

    // QC steps go mostly to preemptable
    withName: 'FASTQC|FASTQ_SCREEN|SAMTOOLS_FLAGSTAT|SAMTOOLS_STATS|SAMTOOLS_IDXSTATS|PICARD_COLLECTMULTIPLEMETRICS|QUALIMAP_BAMQC' {
        // 70% preemptable, 30% cpu
        queue = { 
            new Random().nextInt(100) < 70 ? 'preemptable' : '$SLURM_QUEUE' 
        }
    
	// Alternate on retry
    	clusterOptions = { 
            def partition = task.attempt % 2 == 1 ? 'preemptable' : '$SLURM_QUEUE'
            "--partition=${partition}"
    	}
    
	errorStrategy = 'retry'
    	maxRetries = 10
    }

}

// Executor Configuration
executor {
    \$slurm {
        queueSize = $MAX_JOBS
        submitRateLimit = '$SUBMIT_RATE'
        pollInterval = '30 sec'
    }
}

// Reporting Configuration
report {
    enabled = true
    file = "\${params.outdir}/../logs/reports/nextflow_report.html"
    overwrite = true
}

timeline {
    enabled = true
    file = "\${params.outdir}/../logs/reports/timeline_report.html"
    overwrite = true
}

dag {
    enabled = true
    file = "\${params.outdir}/../logs/reports/pipeline_dag.svg"
    overwrite = true
}

trace {
    enabled = true
    file = "\${params.outdir}/../logs/reports/trace.txt"
    overwrite = true
}

// Manifest
manifest {
    name = '$PROJECT_NAME'
    description = '$PROJECT_DESCRIPTION'
    version = '1.0.0'
    nextflowVersion = '>=23.04.0'
    defaultBranch = 'main'
}
EOL

echo -e "${GREEN}✓ Generated Nextflow config at $CONFIG_DIR/nextflow.config${NC}"

# Generate samples.csv
echo -e "\n${BLUE}Generating sample sheet...${NC}"

if find_sample_ids; then
    # Create samples.csv with detected samples
    cat > "$CONFIG_DIR/samples.csv" << EOL
sample_id,description,group,batch,options
EOL
    
    SAMPLE_COUNT=1
    while read -r sample_id; do
        if [ -n "$sample_id" ]; then
            echo "$sample_id,${sample_id}_description,group1,batch1," >> "$CONFIG_DIR/samples.csv"
            ((SAMPLE_COUNT++))
        fi
    done <<< "$SAMPLE_IDS"
    
    echo -e "${GREEN}✓ Generated sample sheet with $(($SAMPLE_COUNT-1)) samples at $CONFIG_DIR/samples.csv${NC}"
else
    # Create samples.csv with example entries
    cat > "$CONFIG_DIR/samples.csv" << EOL
sample_id,description,group,batch,options
sample1,Sample 1 description,control,batch1,
sample2,Sample 2 description,treatment,batch1,
sample3,Sample 3 description,control,batch2,
sample4,Sample 4 description,treatment,batch2,
EOL
    echo -e "${YELLOW}✓ Generated example sample sheet at $CONFIG_DIR/samples.csv${NC}"
    echo -e "${YELLOW}  Please edit this file to match your actual samples${NC}"
fi

# Create run_pipeline.sh
cat > "$PROJECT_DIR/run_pipeline.sh" << 'EOL'
#!/bin/bash

# ================================
# Pipeline Execution Script
# ================================
# Project: PROJECT_NAME_PLACEHOLDER
# Generated: DATE_PLACEHOLDER
# ================================

#SBATCH --job-name=PROJECT_NAME_PLACEHOLDER
#SBATCH --output=logs/pipeline-%j.log
#SBATCH --error=logs/pipeline-%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --account=SLURM_ACCOUNT_PLACEHOLDER

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'


# Load required modules
echo -e "${YELLOW}Loading modules...${NC}"
if [[ "${CONDA_DEFAULT_ENV:-}" != "single-cell-cnv" ]]; then
    echo "Activating single-cell-cnv environment..."
    conda activate single-cell-cnv || {
        echo "Error: Failed to activate single-cell-cnv environment"
        exit 1
    }
fi

# Check if config files exist
if [ ! -f "config/nextflow.config" ]; then
    echo -e "${RED}Error: config/nextflow.config not found${NC}"
    echo -e "${YELLOW}Please run: nf_cnv/init.sh${NC}"
    exit 1
fi

# Run the pipeline
nextflow run nf_cnv/main.nf -c config/nextflow.config -resume "$@"


# Check exit status
EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
    echo -e "${RED}Pipeline failed with exit code: $EXIT_CODE${NC}"
    echo -e "${YELLOW}Check logs in: logs/${NC}"
fi

echo -e "${BLUE}Finished at: $(date)${NC}"
exit $EXIT_CODE
EOL

# Replace placeholders in run_pipeline.sh
sed -i "s/PROJECT_NAME_PLACEHOLDER/$PROJECT_NAME/g" "$PROJECT_DIR/run_pipeline.sh"
sed -i "s/DATE_PLACEHOLDER/$(date)/g" "$PROJECT_DIR/run_pipeline.sh"
sed -i "s/SLURM_ACCOUNT_PLACEHOLDER/$SLURM_ACCOUNT/g" "$PROJECT_DIR/run_pipeline.sh"

chmod +x "$PROJECT_DIR/run_pipeline.sh"
echo -e "${GREEN}✓ Created executable run script at $PROJECT_DIR/run_pipeline.sh${NC}"

# Create comprehensive README.md
cat > "$PROJECT_DIR/README.md" << EOL
# $PROJECT_NAME - Single Cell CNV Analysis

## Project Overview

**Project Name:** $PROJECT_NAME  
**Description:** $PROJECT_DESCRIPTION  
**Author:** $PROJECT_AUTHOR  
**Contact:** $PROJECT_EMAIL  
**Created:** $(date)  
**Pipeline:** nf-cnv (Single Cell Copy Number Variation Analysis)

## Directory Structure

\`\`\`
$PROJECT_NAME/
├── config/                    # Configuration files
│   ├── nextflow.config        # Main pipeline configuration
│   └─ samples.csv             # Sample metadata
├── seqdata/                   # Raw sequencing data (FASTQ files)
├── results/                   # Pipeline outputs
│   ├── merged/                # Concatenated FASTQ files
│   ├── wsplit/                # Barcode-split cell files
│   ├── bam_out_md/            # Bam File with marked duplicates
│   ├── bam_out_dd/            # Deduplicated bam file
│   ├── bam_out_{fw|rv}        # Mate-Clipped Bam file fw:forward read, rv:reverse read
│   ├── wsplit/                # Barcode-split cell files
│   └── qc/                    # QC statistics and reports
├── logs/                      # Execution logs and reports
│   ├── reports/               # Nextflow HTML reports
│   └── pipeline-*.log         # SLURM job logs
├── docs/                      # Documentation
├── nf_cnv/                    # Pipeline source code
├── run_pipeline.sh            # Main execution script
└── README.md                  # This file
\`\`\`

## Pipeline Workflow

The nf-cnv pipeline performs the following steps:

1. **FASTQ Concatenation**: Merges FASTQ files from multiple sequencing runs/lanes per sample
2. **Barcode Splitting**: Splits reads by cell barcode using optimized barcode matching
3. **Quality Control**: Generates comprehensive QC metrics and filtering
4. **Statistics Generation**: Produces barcode distribution and read count statistics
5. **Report Generation**: Creates HTML reports and visualization

## Quick Start

### 1. Review Configuration

Edit the sample sheet to match your data:
\`\`\`bash
emacs -nw config/samples.csv
\`\`\`

Review pipeline parameters:
\`\`\`bash
emacs -nw config/nextflow.config
\`\`\`

### 2. Run the Pipeline

**Interactive execution:**
\`\`\`bash
./run_pipeline.sh
\`\`\`

**Submit to SLURM:**
\`\`\`bash
sbatch run_pipeline.sh
\`\`\`

**Resume a failed run:**
\`\`\`bash
./run_pipeline.sh --resume
\`\`\`

**Dry run (test without execution):**
\`\`\`bash
./run_pipeline.sh --dry-run
\`\`\`

## Configuration

### Sample Sheet Format

The \`config/samples.csv\` file should contain:

| Column | Description | Example |
|--------|-------------|---------|
| sample_id | Unique sample identifier | sample1 |
| description | Sample description | Control_replicate_1 |
| group | Experimental group | control |
| batch | Sequencing batch | batch1 |
| options | Additional options | (optional) |

### Key Parameters

Edit \`config/nextflow.config\` to modify:

- **Barcode settings**: Mismatch tolerance, search length
- **Resource allocation**: CPU, memory per process
- **Quality filters**: Minimum reads per cell
- **Output options**: File formats, intermediate files

## Results

### Output Files

- **\`results/merged/\`**: Sample-level concatenated FASTQ files
- **\`results/wsplit/\`**: Cell-specific FASTQ files after barcode demultiplexing
- **\`results/stats/\`**: QC metrics, barcode statistics, and summary reports

### Reports

- **\`logs/reports/nextflow_report.html\`**: Execution summary and resource usage
- **\`logs/reports/timeline_report.html\`**: Timeline of process execution
- **\`logs/reports/pipeline_dag.svg\`**: Pipeline workflow diagram

## Troubleshooting

### Common Issues

1. **Pipeline fails to start**
   - Check that \`config/samples.csv\` exists and is properly formatted
   - Verify FASTQ files are present in \`seqdata/\` directory

2. **SLURM job fails**
   - Check SLURM account and partition settings in \`config/nextflow.config\`
   - Review resource requirements (memory, CPU)

3. **Barcode splitting issues**
   - Verify barcode file path and format
   - Adjust mismatch tolerance and search parameters

### Getting Help

- Check execution logs in \`logs/\` directory
- Review Nextflow reports for resource usage and errors
- Examine process-specific logs for detailed error messages

## Advanced Usage

### Custom Profiles

Create custom execution profiles for different environments:

\`\`\`bash
./run_pipeline.sh --profile custom_profile
\`\`\`

### Resource Optimization

For large datasets, consider adjusting:
- Increase memory allocation for barcode splitting
- Adjust concurrent job limits based on cluster capacity
- Enable intermediate file cleanup to save storage

### Integration

Pipeline outputs can be integrated with:
- Single cell analysis workflows (Seurat, Scanpy)
- CNV detection tools
- Downstream statistical analysis

## Project Notes

### TODO
- [ ] Validate sample sheet entries
- [ ] Optimize barcode parameters for your data
- [ ] Review quality control thresholds
- [ ] Set up downstream analysis workflows

### Changes Log
- $(date): Project initialized with nf-cnv pipeline

---

For questions about this project, contact: $PROJECT_EMAIL
EOL

echo -e "${GREEN}✓ Created comprehensive README at $PROJECT_DIR/README.md${NC}"

# Create a project configuration summary
cat > "$CONFIG_DIR/project_info.yaml" << EOL
# Project Configuration Summary
# Generated: $(date)

project:
  name: "$PROJECT_NAME"
  description: "$PROJECT_DESCRIPTION"
  author: "$PROJECT_AUTHOR"
  email: "$PROJECT_EMAIL"
  created: "$(date)"

directories:
  seqdata: "$SEQDATA_DIR"
  config: "$CONFIG_DIR"
  results: "$RESULTS_DIR"
  logs: "$LOGS_DIR"

barcodes:
  file: "$BARCODE_FILE"
  max_mismatches: $BARCODE_MISMATCHES
  search_length: $BARCODE_SEARCH_LEN
  reverse_complement: $BARCODE_RC
  trim_barcode: $BARCODE_TRIM
EOL

echo -e "${GREEN}✓ Created project info at $CONFIG_DIR/project_info.yaml${NC}"

# Create a quick reference guide
cat > "$PROJECT_DIR/QUICK_REFERENCE.md" << EOL
# Quick Reference - $PROJECT_NAME

## Essential Commands

\`\`\`bash
# Run pipeline
./run_pipeline.sh

# Submit to SLURM
sbatch run_pipeline.sh

# Resume failed run
./run_pipeline.sh --resume

# Test run (no execution)
./run_pipeline.sh --dry-run

# Check job status
squeue -u \$USER

# Cancel job
scancel <job_id>
\`\`\`

## Important Files

- **\`config/samples.csv\`** - Edit your sample information here
- **\`config/nextflow.config\`** - Pipeline parameters and resources
- **\`logs/pipeline-*.log\`** - Check for errors here
- **\`results/\`** - Your output files will be here

## Quick Checks

\`\`\`bash
# Verify sample sheet
head -5 config/samples.csv

# Check for FASTQ files
find seqdata/ -name "*.fastq.gz" | head -5

# Monitor disk usage
du -sh results/ logs/

# View latest log
tail -f logs/pipeline-*.log
\`\`\`

## Pipeline Status

\`\`\`bash
# Check if pipeline is running
ps aux | grep nextflow

# View resource usage
sstat <job_id>

# Check node information
sinfo -N -l
\`\`\`
EOL

echo -e "${GREEN}✓ Created quick reference at $PROJECT_DIR/QUICK_REFERENCE.md${NC}"

# Create .gitignore
cat > "$PROJECT_DIR/.gitignore" << EOL
# Results and intermediate files
results/
logs/
.nextflow*
work/

# Temporary files
*.tmp
*.temp
*~

# Data files (comment out if you want to track them)
seqdata/
*.fastq.gz
*.fq.gz

# System files
.DS_Store
Thumbs.db

# Editor files
*.swp
*.swo
*~
EOL

echo -e "${GREEN}✓ Created .gitignore file${NC}"

# Final summary and instructions
echo -e "\n${CYAN}=================================================${NC}"
echo -e "${GREEN} Project '$PROJECT_NAME' Successfully Initialized! ${NC}"
echo -e "${CYAN}=================================================${NC}"

echo -e "\n${BLUE} Project Structure Created:${NC}"
echo -e "   ${GREEN}✓${NC} Configuration files in config/"
echo -e "   ${GREEN}✓${NC} Execution script: run_pipeline.sh"
echo -e "   ${GREEN}✓${NC} Documentation: README.md"
echo -e "   ${GREEN}✓${NC} Quick reference guide"
echo -e "   ${GREEN}✓${NC} Directory structure for results and logs"

echo -e "\n${PURPLE} Next Steps:${NC}"
echo -e "   ${YELLOW}1.${NC} Review and edit sample information:"
echo -e "      ${CYAN}nano config/samples.csv${NC}"
echo -e "   ${YELLOW}2.${NC} Verify pipeline configuration:"
echo -e "      ${CYAN}nano config/nextflow.config${NC}"
echo -e "   ${YELLOW}3.${NC} Place your FASTQ files in:"
echo -e "      ${CYAN}$SEQDATA_DIR${NC}"
echo -e "   ${YELLOW}4.${NC} Test the pipeline:"
echo -e "      ${CYAN}./run_pipeline.sh --dry-run${NC}"
echo -e "   ${YELLOW}5.${NC} Run the pipeline:"
echo -e "      ${CYAN}sbatch run_pipeline.sh${NC}"

echo -e "\n${BLUE} Documentation:${NC}"
echo -e "   • ${CYAN}README.md${NC} - Complete project documentation"
echo -e "   • ${CYAN}QUICK_REFERENCE.md${NC} - Essential commands"
echo -e "   • ${CYAN}config/project_info.yaml${NC} - Configuration summary"

echo -e "\n${GREEN}Happy single-cell analysis! 🧬✨${NC}"
echo -e "${CYAN}=================================================${NC}"
