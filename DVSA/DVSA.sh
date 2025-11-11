#!/bin/bash

# ===================================================================================
# Differential viral segment abundence analysis (DVSA) pipeline for RNA- or sRNA-seq
# Methods of ViralScout methodology
# Users can define any their pre-assembled contigs as viral benchmarks by providing names
# Caculate the Euclidean distance (TRM, FC) between unknown contigs with the centroid of viral benchmarks
# Please put all needed files in current directory
# ===================================================================================


set -e
set -u
set -o pipefail

# ---------- Logging utilities ----------
log_info()      { echo "[INFO] $(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_error()     { echo "[ERROR] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2; }
log_warning()   { echo "[WARNING] $(date '+%Y-%m-%d %H:%M:%S') - $1" >&2; }
log_success()   { echo "[SUCCESS] $(date '+%Y-%m-%d %H:%M:%S') - $1"; }
log_step() {
    echo ""
    echo "===================================================================="
    echo "[STEP] $(date '+%Y-%m-%d %H:%M:%S') - $1"
    echo "===================================================================="
}

# ---------- Command runner ----------
run_command() {
    cmd="$1"
    step_name="$2"
    log_info "Executing: $cmd"
    if eval "$cmd"; then
        log_success "$step_name completed successfully"
    else
        exit_code=$?
        log_error "$step_name failed with exit code $exit_code"
        log_error "Failed command: $cmd"
        exit $exit_code
    fi
}

# ---------- File checker ----------
check_file() {
    file="$1"
    desc="$2"
    if [[ ! -f "$file" ]]; then
        log_error "$desc file does not exist: $file"
        return 1
    fi
    file_size=$(stat -c%s "$file" 2>/dev/null || stat -f%z "$file" 2>/dev/null)
    if [[ $file_size -eq 0 ]]; then
        log_warning "$desc file is empty: $file"
    fi
    log_info "$desc: $file (size: ${file_size} bytes)"
    return 0
}

# ---------- Detect read type ----------
detect_read_type() {
    local positive_file="$1"
    local negative_file="$2"
    
    # Detect files
    local positive_exists=true
    local negative_exists=true
    
    # Detect virus-infected file
    if [[ "$positive_file" == *","* ]]; then
        IFS=',' read -ra pos_files <<< "$positive_file"
        for file in "${pos_files[@]}"; do
            if [[ ! -f "$file" ]]; then
                log_error "Positive sample file does not exist: $file"
                positive_exists=false
            fi
        done
    else
        if [[ ! -f "$positive_file" ]]; then
            log_error "Positive sample file does not exist: $positive_file"
            positive_exists=false
        fi
    fi
    
    # Detect virus-uninfected file
    if [[ "$negative_file" == *","* ]]; then
        IFS=',' read -ra neg_files <<< "$negative_file"
        for file in "${neg_files[@]}"; do
            if [[ ! -f "$file" ]]; then
                log_error "Negative sample file does not exist: $file"
                negative_exists=false
            fi
        done
    else
        if [[ ! -f "$negative_file" ]]; then
            log_error "Negative sample file does not exist: $negative_file"
            negative_exists=false
        fi
    fi
    
    if [[ "$positive_exists" == false ]] || [[ "$negative_exists" == false ]]; then
        return 1
    fi
    
    # Detect paired-end or not
    if [[ "$positive_file" == *"_R1."* ]] && [[ "$positive_file" == *"_R2."* ]]; then
        echo "paired"
    elif [[ "$negative_file" == *"_R1."* ]] && [[ "$negative_file" == *"_R2."* ]]; then
        echo "paired"
    elif [[ "$positive_file" == *"R1"* ]] && [[ -f "${positive_file/R1/R2}" ]]; then
        echo "paired"
    elif [[ "$negative_file" == *"R1"* ]] && [[ -f "${negative_file/R1/R2}" ]]; then
        echo "paired"
    else
        echo "single"
    fi
}

# ---------- Help message ----------
usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required parameters:"
    echo "  --contigs FILE          Reference contigs fasta file"
    echo "  --positive FILE         Positive sample fastq file(s)"
    echo "  --negative FILE         Negative sample fastq file(s)"
    echo "  --benchmark FILE        Benchmark contigs file"
    echo "  --output DIR            Output directory"
    echo ""
    echo "For paired-end reads, use:"
    echo "  --positive FILE1,FILE2  Comma-separated R1 and R2 files"
    echo "  --negative FILE1,FILE2  Comma-separated R1 and R2 files"
    echo ""
    echo "Performance options:"
    echo "  --threads INT           Number of threads (default: 32)"
    echo "  --bowtie2-threads INT   Threads for bowtie2 (default: 32)"
    echo "  --samtools-threads INT  Threads for samtools (default: 32)"
    echo ""
    echo "Analysis options:"
    echo "  --nearest-percent FLOAT Percentage of nearest contigs to select (default: 0.1)"
    echo ""
    echo "Other options:"
    echo "  --force                 Force rerun from beginning"
    echo "  -h, --help              Show help message"
}

# ---------- Default parameters ----------
CONTIGS_FILE=""
POSITIVE_READS=""
NEGATIVE_READS=""
BENCHMARK_FILE=""
OUTPUT_DIR=""
THREADS=32
BOWTIE2_THREADS=32
SAMTOOLS_THREADS=32
NEAREST_PERCENT=0.1
FORCE_RERUN=false

# ---------- Parse arguments ----------
while [[ $# -gt 0 ]]; do
    case $1 in
        # Required parameters
        --contigs) CONTIGS_FILE="$2"; shift 2 ;;
        --positive) POSITIVE_READS="$2"; shift 2 ;;
        --negative) NEGATIVE_READS="$2"; shift 2 ;;
        --benchmark) BENCHMARK_FILE="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        
        # Performance options
        --threads) THREADS="$2"; shift 2 ;;
        --bowtie2-threads) BOWTIE2_THREADS="$2"; shift 2 ;;
        --samtools-threads) SAMTOOLS_THREADS="$2"; shift 2 ;;
        
        # Analysis options
        --nearest-percent) NEAREST_PERCENT="$2"; shift 2 ;;
        
        # Other options
        --force) FORCE_RERUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

# ---------- Check required inputs ----------
if [[ -z "$CONTIGS_FILE" || -z "$POSITIVE_READS" || -z "$NEGATIVE_READS" || -z "$BENCHMARK_FILE" || -z "$OUTPUT_DIR" ]]; then
    log_error "Missing required parameters"
    usage
    exit 1
fi

# ---------- Convert relative paths to absolute paths ----------
convert_to_absolute_path() {
    local path="$1"
    if [[ "$path" == /* ]]; then
        echo "$path"
    else
        echo "$(pwd)/$path"
    fi
}

CONTIGS_FILE=$(convert_to_absolute_path "$CONTIGS_FILE")
BENCHMARK_FILE=$(convert_to_absolute_path "$BENCHMARK_FILE")

# Process FASTQ file
convert_comma_separated_paths() {
    local paths="$1"
    if [[ "$paths" == *","* ]]; then
        IFS=',' read -ra files <<< "$paths"
        local result=""
        for file in "${files[@]}"; do
            if [[ -n "$result" ]]; then
                result="$result,"
            fi
            result="$result$(convert_to_absolute_path "$file")"
        done
        echo "$result"
    else
        convert_to_absolute_path "$paths"
    fi
}

POSITIVE_READS=$(convert_comma_separated_paths "$POSITIVE_READS")
NEGATIVE_READS=$(convert_comma_separated_paths "$NEGATIVE_READS")

# ---------- Check if input files exist ----------
log_step "Checking input files"
check_file "$CONTIGS_FILE" "Contigs file"
check_file "$BENCHMARK_FILE" "Benchmark file"

# Detect FASTQ file
if [[ "$POSITIVE_READS" == *","* ]]; then
    IFS=',' read -ra POSITIVE_FILES <<< "$POSITIVE_READS"
    for file in "${POSITIVE_FILES[@]}"; do
        check_file "$file" "Positive sample file"
    done
else
    check_file "$POSITIVE_READS" "Positive sample file"
fi

if [[ "$NEGATIVE_READS" == *","* ]]; then
    IFS=',' read -ra NEGATIVE_FILES <<< "$NEGATIVE_READS"
    for file in "${NEGATIVE_FILES[@]}"; do
        check_file "$file" "Negative sample file"
    done
else
    check_file "$NEGATIVE_READS" "Negative sample file"
fi

# ---------- Setup output directory ----------
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit 1

LOG_FILE="analysis_pipeline_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -i "$LOG_FILE")
exec 2>&1

log_step "Starting Seq Analysis Pipeline"
log_info "Contigs file: $CONTIGS_FILE"
log_info "Positive reads: $POSITIVE_READS"
log_info "Negative reads: $NEGATIVE_READS"
log_info "Benchmark file: $BENCHMARK_FILE"
log_info "Output directory: $OUTPUT_DIR"
log_info "Threads: $THREADS"
log_info "Nearest contigs percentage: ${NEAREST_PERCENT}%"
log_info "Force rerun: $FORCE_RERUN"

# ---------- Detect read type ----------
log_step "Detecting read type"
READ_TYPE=$(detect_read_type "$POSITIVE_READS" "$NEGATIVE_READS")
if [[ $? -ne 0 ]]; then
    log_error "Failed to detect read type due to missing files"
    exit 1
fi
log_info "Detected read type: $READ_TYPE"

# Parse file names for paired-end reads
if [[ "$READ_TYPE" == "paired" ]]; then
    IFS=',' read -ra POSITIVE_FILES <<< "$POSITIVE_READS"
    IFS=',' read -ra NEGATIVE_FILES <<< "$NEGATIVE_READS"
    
    if [[ ${#POSITIVE_FILES[@]} -ne 2 ]] || [[ ${#NEGATIVE_FILES[@]} -ne 2 ]]; then
        log_error "Paired-end mode requires exactly 2 files for both positive and negative samples"
        exit 1
    fi
    
    POSITIVE_R1="${POSITIVE_FILES[0]}"
    POSITIVE_R2="${POSITIVE_FILES[1]}"
    NEGATIVE_R1="${NEGATIVE_FILES[0]}"
    NEGATIVE_R2="${NEGATIVE_FILES[1]}"
    
    log_info "Positive R1: $POSITIVE_R1"
    log_info "Positive R2: $POSITIVE_R2"
    log_info "Negative R1: $NEGATIVE_R1"
    log_info "Negative R2: $NEGATIVE_R2"
else
    log_info "Positive file: $POSITIVE_READS"
    log_info "Negative file: $NEGATIVE_READS"
fi

# ---------- Define output file names ----------
INDEX_PREFIX="contigs_index"
POSITIVE_SAM="positive_sample.sam"
POSITIVE_BAM="positive_sample.bam"
POSITIVE_SORTED_BAM="positive_sample_sorted.bam"
NEGATIVE_SAM="negative_sample.sam"
NEGATIVE_BAM="negative_sample.bam"
NEGATIVE_SORTED_BAM="negative_sample_sorted.bam"
CONTIG_LENGTHS="contig_lengths.txt"

# ---------- Check if we can resume from sorted BAM files ----------
if [[ "$FORCE_RERUN" = false ]] && [[ -f "$POSITIVE_SORTED_BAM" ]] && [[ -f "$NEGATIVE_SORTED_BAM" ]] && [[ -f "${POSITIVE_SORTED_BAM}.bai" ]] && [[ -f "${NEGATIVE_SORTED_BAM}.bai" ]]; then
    log_success "Found existing sorted BAM files, resuming from analysis"
    SKIP_ALIGNMENT=true
else
    log_info "Starting full analysis pipeline from beginning"
    SKIP_ALIGNMENT=false
fi

# ---------- Step 1: Build Bowtie2 index ----------
if [[ "$SKIP_ALIGNMENT" = false ]]; then
    log_step "Step 1: Building Bowtie2 index"
    
    if [[ ! -f "${INDEX_PREFIX}.1.bt2" ]] || [[ "$FORCE_RERUN" = true ]]; then
        run_command "bowtie2-build --threads $BOWTIE2_THREADS '$CONTIGS_FILE' '$INDEX_PREFIX'" "Index building"
        check_file "${INDEX_PREFIX}.1.bt2" "Bowtie2 index"
    else
        log_success "Bowtie2 index already exists, skipping"
    fi
fi

# ---------- Step 2: Alignment ----------
if [[ "$SKIP_ALIGNMENT" = false ]]; then
    log_step "Step 2: Sequence alignment"
    
    # Positive sample alignment
    if [[ "$READ_TYPE" == "paired" ]]; then
        run_command "bowtie2 -x '$INDEX_PREFIX' -1 '$POSITIVE_R1' -2 '$POSITIVE_R2' -S '$POSITIVE_SAM' --very-sensitive -p $BOWTIE2_THREADS --no-unal" "Positive sample alignment (paired-end)"
    else
        run_command "bowtie2 -x '$INDEX_PREFIX' -U '$POSITIVE_READS' -S '$POSITIVE_SAM' -N 0 -L 15 --very-sensitive -p $BOWTIE2_THREADS --no-unal" "Positive sample alignment (single-end)"
    fi
    
    # Negative sample alignment
    if [[ "$READ_TYPE" == "paired" ]]; then
        run_command "bowtie2 -x '$INDEX_PREFIX' -1 '$NEGATIVE_R1' -2 '$NEGATIVE_R2' -S '$NEGATIVE_SAM' --very-sensitive -p $BOWTIE2_THREADS --no-unal" "Negative sample alignment (paired-end)"
    else
        run_command "bowtie2 -x '$INDEX_PREFIX' -U '$NEGATIVE_READS' -S '$NEGATIVE_SAM' -N 0 -L 15 --very-sensitive -p $BOWTIE2_THREADS --no-unal" "Negative sample alignment (single-end)"
    fi
    
    check_file "$POSITIVE_SAM" "Positive sample SAM"
    check_file "$NEGATIVE_SAM" "Negative sample SAM"
fi

# ---------- Step 3: SAM to BAM conversion ----------
if [[ "$SKIP_ALIGNMENT" = false ]]; then
    log_step "Step 3: SAM to BAM conversion"
    
    run_command "samtools view -b --threads $SAMTOOLS_THREADS '$POSITIVE_SAM' > '$POSITIVE_BAM'" "Convert positive SAM to BAM"
    run_command "samtools view -b --threads $SAMTOOLS_THREADS '$NEGATIVE_SAM' > '$NEGATIVE_BAM'" "Convert negative SAM to BAM"
    
    check_file "$POSITIVE_BAM" "Positive sample BAM"
    check_file "$NEGATIVE_BAM" "Negative sample BAM"
fi

# ---------- Step 4: Sort and index BAM files ----------
if [[ "$SKIP_ALIGNMENT" = false ]]; then
    log_step "Step 4: Sort and index BAM files"
    
    run_command "samtools sort --threads $SAMTOOLS_THREADS '$POSITIVE_BAM' -o '$POSITIVE_SORTED_BAM'" "Sort positive BAM"
    run_command "samtools sort --threads $SAMTOOLS_THREADS '$NEGATIVE_BAM' -o '$NEGATIVE_SORTED_BAM'" "Sort negative BAM"
    run_command "samtools index '$POSITIVE_SORTED_BAM'" "Index positive BAM"
    run_command "samtools index '$NEGATIVE_SORTED_BAM'" "Index negative BAM"
    
    check_file "$POSITIVE_SORTED_BAM" "Sorted positive BAM"
    check_file "$NEGATIVE_SORTED_BAM" "Sorted negative BAM"
    check_file "${POSITIVE_SORTED_BAM}.bai" "Positive BAM index"
    check_file "${NEGATIVE_SORTED_BAM}.bai" "Negative BAM index"
    
    # Clean up intermediate files
    log_step "Cleaning intermediate files"
    rm -f "$POSITIVE_SAM" "$NEGATIVE_SAM" "$POSITIVE_BAM" "$NEGATIVE_BAM"
    log_success "Cleaned intermediate SAM and unsorted BAM files"
fi

# ---------- Step 5: Generate contig lengths file ----------
log_step "Step 5: Generating contig lengths file"

if [[ ! -f "$CONTIG_LENGTHS" ]] || [[ "$FORCE_RERUN" = true ]]; then
    run_command "samtools faidx '$CONTIGS_FILE'" "Index contigs file"
    run_command "awk '{print \$1\"\t\"\$2}' '${CONTIGS_FILE}.fai' > '$CONTIG_LENGTHS'" "Extract contig lengths"
    check_file "$CONTIG_LENGTHS" "Contig lengths"
else
    log_success "Contig lengths file already exists, skipping"
fi

# ---------- Step 6: Run seq analysis ----------
log_step "Step 6: Running seq differential expression analysis"

ANALYSIS_OUTPUT="seq_analysis"

run_command "python3 ../DVSAanalyzer.py --positive_bam '$POSITIVE_SORTED_BAM' --negative_bam '$NEGATIVE_SORTED_BAM' --contig_lengths '$CONTIG_LENGTHS' --benchmark_contigs '$BENCHMARK_FILE' --output '$ANALYSIS_OUTPUT' --nearest_percent '$NEAREST_PERCENT'" "seq analysis" 

# ---------- Final summary ----------
log_step "Pipeline completion report"
log_success "Seq Analysis Pipeline Completed Successfully!"
echo ""
echo "Output files:"
echo "  - Alignment files:"
echo "    - $POSITIVE_SORTED_BAM"
echo "    - $NEGATIVE_SORTED_BAM"
echo "  - Analysis results:"
echo "    - ${ANALYSIS_OUTPUT}_seq_data.csv"
echo "    - ${ANALYSIS_OUTPUT}_with_distances.csv"
echo "    - ${ANALYSIS_OUTPUT}_nearest_contigs.csv"
echo "    - ${ANALYSIS_OUTPUT}_feature_space.svg"
echo "  - Log file: $LOG_FILE"
echo ""
echo "To rerun the analysis, use:"
echo "  $0 --contigs '$CONTIGS_FILE' --positive '$POSITIVE_READS' --negative '$NEGATIVE_READS' --benchmark '$BENCHMARK_FILE' --output '$OUTPUT_DIR' --threads $THREADS --nearest-percent $NEAREST_PERCENT --force"
