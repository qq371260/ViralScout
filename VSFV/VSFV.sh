#!/bin/bash

# ==================================================================================================
# Viral small-RNA feature vector (VSFV) analysis - MODULAR VERSION
# A method of ViralScout methodology
# It uses variable-dimensional feature vectors counted from uniquely mapped sRNAs of each contig
# Feature modes: size (length only), size_P_5nt (length + 5' nt), sizeXstr (length x strand), sizeX5nt (length x 5' nt), sizeX5ntXstr (length x 5' nt x strand)
# Users can define any their pre-assembled contigs as viral benchmarks by providing names
# Other contigs were compared with the benchmarks to find virus-like contigs
# The comparison is based on pairwise correlation coefficients of feature vectors (Spearman or Pearson)
# Please put all needed files in current directory
# ==================================================================================================


set -uo pipefail

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

# ---------- Extract base filename ----------
get_basename() {
    local filename="$1"
    # Remove directory path and file extension
    basename=$(basename "$filename")
    # Remove common fastq/fasta extensions
    basename="${basename%.fastq}"
    basename="${basename%.fq}"
    basename="${basename%.fasta}"
    basename="${basename%.fa}"
    basename="${basename%.fna}"
    echo "$basename"
}

# ---------- Help message ----------
usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Required parameters:"
    echo "  -f, --fasta FILE         Reference genome fasta file (required)"
    echo "  -r, --reads FILE         sRNA sequencing fastq file (required)"
    echo ""
    echo "Input/output options:"
    echo "  -b, --benchmark FILE     Benchmark contigs file (optional)"
    echo "  -o, --output DIR         Output directory (default: ./sRNA_analysis)"
    echo ""
    echo "Performance options:"
    echo "  -t, --threads INT        Total threads (default: all cores)"
    echo "  --bowtie2-threads INT    Threads for bowtie2 (default: 8)"
    echo "  --samtools-threads INT   Threads for samtools (default: 4)"
    echo "  --sort-memory STR        Memory for samtools sort (default: 2G)"
    echo ""
    echo "sRNA analysis options:"
    echo "  --min-length INT         Minimum sRNA length (default: 18)"
    echo "  --max-length INT         Maximum sRNA length (default: 30)"
    echo "  --keep-non-unique        Keep non-uniquely mapped reads (default: filter to unique only)"
    echo ""
    echo "contig_filter.py options:"
    echo "  -n, --number_21_nt INT   21nt threshold for filtering (default: 10)"
    echo ""
    echo "vector.py options:"
    echo "  --vector-cores INT       Cores for vector.py (default: same as samtools-threads)"
    echo "  -m, --mode STR           Feature vector mode: size/size_P_5nt/sizeXstr/sizeX5nt/sizeX5ntXstr (default: sizeX5nt)"
    echo "                           size: length only"
    echo "                           size_P_5nt: length + 5' nucleotides" 
    echo "                           sizeXstr: length x strands"
    echo "                           sizeX5nt: length x 5' nucleotides"
    echo "                           sizeX5ntXstr: length x 5' nucleotides x strands"
    echo ""
    echo "corr.py options:"
    echo "  --corr-cores INT         Cores for corr.py (default: same as samtools-threads)"
    echo "  --pearson                Use Pearson correlation instead of Spearman (default: Spearman)"
    echo "  --mean-r FLOAT          Composite score threshold (default: 0.8)"
    echo "  --p-value FLOAT          P-value threshold (default: 0.05)"
    echo ""
    echo "contig_feature.py options:"
    echo "  --feature-cores INT      Cores for contig_feature.py (default: same as samtools-threads)"
    echo ""
    echo "Other options:"
    echo "  --force                  Force rerun from beginning (ignore existing files)"
    echo "  -h, --help               Show help message"
}

# ---------- Default parameters ----------
FASTA_FILE=""
READS_FILE=""
BENCHMARK_FILE=""
OUTPUT_DIR="./sRNA_analysis"
TOTAL_THREADS=$(nproc)
BOWTIE2_THREADS=8
SAMTOOLS_THREADS=4
SORT_MEMORY="2G"
MIN_LENGTH=18
MAX_LENGTH=30
KEEP_NON_UNIQUE=false
FORCE_RERUN=false
USE_PEARSON=false

# contig_filter.py parameters
NUMBER_21_NT=10

# vector.py parameters  
VECTOR_CORES="$SAMTOOLS_THREADS"

# corr.py parameters
CORR_CORES="$SAMTOOLS_THREADS"
MEAN_R=0.8
P_VALUE=0.05
VECTOR_MODE="sizeX5nt"

# contig_feature.py parameters
FEATURE_CORES="$SAMTOOLS_THREADS"

# ---------- Parse arguments ----------
while [[ $# -gt 0 ]]; do
    case $1 in
        # Required parameters
        -f|--fasta) FASTA_FILE="$2"; shift 2 ;;
        -r|--reads) READS_FILE="$2"; shift 2 ;;
        
        # Input/output options
        -b|--benchmark) BENCHMARK_FILE="$2"; shift 2 ;;
        -o|--output) OUTPUT_DIR="$2"; shift 2 ;;
        
        # Performance options
        -t|--threads) TOTAL_THREADS="$2"; shift 2 ;;
        --bowtie2-threads) BOWTIE2_THREADS="$2"; shift 2 ;;
        --samtools-threads) SAMTOOLS_THREADS="$2"; shift 2 ;;
        --sort-memory) SORT_MEMORY="$2"; shift 2 ;;
        
        # sRNA analysis options
        --min-length) MIN_LENGTH="$2"; shift 2 ;;
        --max-length) MAX_LENGTH="$2"; shift 2 ;;
        --keep-non-unique) KEEP_NON_UNIQUE=true; shift ;;
        
        # contig_filter.py options
        -n|--number_21_nt) NUMBER_21_NT="$2"; shift 2 ;;
        
        # vector.py options
        --vector-cores) VECTOR_CORES="$2"; shift 2 ;;
        -m|--mode) VECTOR_MODE="$2"; shift 2 ;;
        
        # corr.py options
        --corr-cores) CORR_CORES="$2"; shift 2 ;;
        --pearson) USE_PEARSON=true; shift ;;
        --mean-r) MEAN_R="$2"; shift 2 ;;
        --p-value) P_VALUE="$2"; shift 2 ;;
        
        # contig_feature.py options
        --feature-cores) FEATURE_CORES="$2"; shift 2 ;;
        
        # Other options
        --force) FORCE_RERUN=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) log_error "Unknown option: $1"; usage; exit 1 ;;
    esac
done

# ---------- Set dependent defaults ----------
# If vector cores not explicitly set, use samtools threads
if [[ "$VECTOR_CORES" == "$SAMTOOLS_THREADS" ]]; then
    VECTOR_CORES="$SAMTOOLS_THREADS"
fi

# If corr cores not explicitly set, use samtools threads
if [[ "$CORR_CORES" == "$SAMTOOLS_THREADS" ]]; then
    CORR_CORES="$SAMTOOLS_THREADS"
fi

# If feature cores not explicitly set, use samtools threads
if [[ "$FEATURE_CORES" == "$SAMTOOLS_THREADS" ]]; then
    FEATURE_CORES="$SAMTOOLS_THREADS"
fi

# ---------- Check required inputs ----------
if [[ -z "$FASTA_FILE" || -z "$READS_FILE" ]]; then
    log_error "Must provide both fasta and reads files."
    usage
    exit 1
fi

# ---------- Validate parameters ----------
# Validate dimension parameter
if [[ ! "$VECTOR_MODE" =~ ^(size|size_P_5nt|sizeXstr|sizeX5nt|sizeX5ntXstr)$ ]]; then
    log_error "Invalid mode: $VECTOR_MODE. Must be size, size_P_5nt, sizeXstr, sizeX5nt, or sizeX5ntXstr"
    exit 1
fi

# ---------- Set correlation method ----------
if [[ "$USE_PEARSON" = true ]]; then
    CORRELATION_METHOD="pearson"
    CORRELATION_DESCRIPTION="Pearson correlation"
else
    CORRELATION_METHOD="spearman"
    CORRELATION_DESCRIPTION="Spearman correlation"
fi

# ---------- Extract base names for file prefixes ----------
READS_BASENAME=$(get_basename "$READS_FILE")
FASTA_BASENAME=$(get_basename "$FASTA_FILE")
SAMPLE_PREFIX="${READS_BASENAME}_vs_${FASTA_BASENAME}"

# ---------- Setup log file ----------
mkdir -p "$OUTPUT_DIR"
LOG_FILE="${OUTPUT_DIR}/${SAMPLE_PREFIX}_pipeline_$(date '+%Y%m%d_%H%M%S').log"
exec > >(tee -a "$LOG_FILE") 2>&1
log_info "Logging started. All messages will be saved to: $LOG_FILE"

log_info "Sample prefix: $SAMPLE_PREFIX"
log_info "Reads basename: $READS_BASENAME"
log_info "Reference basename: $FASTA_BASENAME"
log_info "Feature vector mode: $VECTOR_MODE"
log_info "Correlation method: $CORRELATION_DESCRIPTION"
log_info "Keep non-unique reads: $KEEP_NON_UNIQUE"

# ---------- Check all required files exist in current directory ----------
log_step "Checking required files in current directory"

# Check benchmark file (if provided) - search in multiple locations
if [[ -n "$BENCHMARK_FILE" ]]; then
    if [[ -f "../$BENCHMARK_FILE" ]]; then
        log_info "Benchmark file found in parent directory: ../$BENCHMARK_FILE"
    elif [[ -f "$BENCHMARK_FILE" ]]; then
        log_info "Benchmark file found in current directory: $BENCHMARK_FILE"
    else
        log_error "Benchmark file not found: $BENCHMARK_FILE"
        log_error "Please check the file exists in either:"
        log_error "  - Current directory: ./$BENCHMARK_FILE"
        log_error "  - Parent directory: ../$BENCHMARK_FILE"
        exit 1
    fi
fi

check_file "$FASTA_FILE" "Reference fasta"
check_file "$READS_FILE" "sRNA reads"

# Check new modular scripts
check_file "contig_filter.py" "Contig filter script"
check_file "vector.py" "Feature vector script"
check_file "corr.py" "Correlation analysis script"
check_file "contig_feature.py" "Contig feature script"

# ---------- Setup output ----------
LOG_FILE="${OUTPUT_DIR}/${SAMPLE_PREFIX}_pipeline_$(date '+%Y%m%d_%H%M%S').log"
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit 1

log_step "Starting MODULAR sRNA Sequencing Data Analysis Pipeline"
log_info "Reference genome: ../$FASTA_FILE"
log_info "Sequencing data: ../$READS_FILE"
log_info "Benchmark file: ${BENCHMARK_FILE:+../$BENCHMARK_FILE}"
log_info "Output directory: $OUTPUT_DIR"
log_info "Sample prefix: $SAMPLE_PREFIX"
log_info "Feature vector mode: $VECTOR_MODE"
log_info "Correlation method: $CORRELATION_DESCRIPTION"
log_info "Thread configuration: Total=$TOTAL_THREADS, bowtie2=$BOWTIE2_THREADS, samtools=$SAMTOOLS_THREADS"
log_info "Sort memory: $SORT_MEMORY"
log_info "sRNA length range: ${MIN_LENGTH}-${MAX_LENGTH} nt"
log_info "21-nt threshold: $NUMBER_21_NT"
log_info "Keep non-unique reads: $KEEP_NON_UNIQUE"
log_info "Force rerun: $FORCE_RERUN"
log_info "All Python scripts confirmed in parent directory"
log_info "Log file: $LOG_FILE"

# ---------- Define other output files ----------
FEATURE_CSV="${SAMPLE_PREFIX}_${VECTOR_MODE}_features.csv"
CORRELATION_PREFIX="${SAMPLE_PREFIX}_${CORRELATION_METHOD}_${VECTOR_MODE}"

# ---------- Check if we can resume from indexed BAM ----------
if [[ "$KEEP_NON_UNIQUE" = true ]]; then
    CHECK_BAM="${SAMPLE_PREFIX}_all_mapped.bam"
    CHECK_BAM_INDEX="${CHECK_BAM}.bai"
else
    CHECK_BAM="${SAMPLE_PREFIX}_unique.bam"
    CHECK_BAM_INDEX="${CHECK_BAM}.bai"
fi

if [[ "$FORCE_RERUN" = false ]] && [[ -f "$CHECK_BAM" ]] && [[ -f "$CHECK_BAM_INDEX" ]]; then
    log_success "Found existing indexed BAM file, resuming from downstream analysis"
    SKIP_TO_DOWNSTREAM=true
    FINAL_BAM="$CHECK_BAM"
    FINAL_BAM_INDEX="$CHECK_BAM_INDEX"
else
    log_info "Starting full analysis pipeline from beginning"
    SKIP_TO_DOWNSTREAM=false
    if [[ "$KEEP_NON_UNIQUE" = true ]]; then
        FINAL_BAM="${SAMPLE_PREFIX}_all_mapped.bam"
        FINAL_BAM_INDEX="${FINAL_BAM}.bai"
    else
        FINAL_BAM="${SAMPLE_PREFIX}_unique.bam"
        FINAL_BAM_INDEX="${FINAL_BAM}.bai"
    fi
    [[ "$FORCE_RERUN" = true ]] && log_info "  - FORCE_RERUN is true"
    [[ ! -f "$CHECK_BAM" ]] && log_info "  - BAM file missing: $CHECK_BAM"
    [[ ! -f "$CHECK_BAM_INDEX" ]] && log_info "  - BAM index missing: $CHECK_BAM_INDEX"
fi

# ---------- Debug information ----------
log_info "Using BAM file: $FINAL_BAM"
log_info "BAM index: $FINAL_BAM_INDEX"
log_info "Skip to downstream: $SKIP_TO_DOWNSTREAM"

if [[ "$SKIP_TO_DOWNSTREAM" = false ]]; then
    # ---------- Step 1: Build Bowtie2 index ----------
    log_step "Step 1: Building Bowtie2 index"
    INDEX_PREFIX="${FASTA_BASENAME}_index"
    run_command "bowtie2-build --threads $BOWTIE2_THREADS '../$FASTA_FILE' '$INDEX_PREFIX'" "Index building"
    check_file "${INDEX_PREFIX}.1.bt2" "Bowtie2 index base file"

    # ---------- Step 2: Alignment ----------
    log_step "Step 2: Sequence alignment"
    BAM_PREFIX="$SAMPLE_PREFIX"
    run_command "bowtie2 --threads $BOWTIE2_THREADS -x '$INDEX_PREFIX' -U '../$READS_FILE' -S '${BAM_PREFIX}.sam' -k 2 -N 0 --no-unal" "Sequence alignment"
    check_file "${BAM_PREFIX}.sam" "Alignment SAM file"

    # ---------- Step 3: SAM to BAM, sorting, indexing ----------
    log_step "Step 3: SAM processing"
    run_command "samtools view -h -F 4 --threads $SAMTOOLS_THREADS '${BAM_PREFIX}.sam' > '${BAM_PREFIX}_mapped.sam'" "Extract mapped reads"
    run_command "samtools view -b --threads $SAMTOOLS_THREADS '${BAM_PREFIX}_mapped.sam' > '${BAM_PREFIX}_mapped.bam'" "Convert SAM to BAM"
    run_command "samtools sort -m $SORT_MEMORY --threads $SAMTOOLS_THREADS '${BAM_PREFIX}_mapped.bam' -o '${BAM_PREFIX}_mapped_sorted.bam'" "Sort BAM file"
    run_command "samtools index -@ $SAMTOOLS_THREADS '${BAM_PREFIX}_mapped_sorted.bam'" "Index BAM file"
    check_file "${BAM_PREFIX}_mapped_sorted.bam" "Sorted BAM file"
    check_file "${BAM_PREFIX}_mapped_sorted.bam.bai" "BAM index file"

    # Clean up intermediates
    log_step "Cleaning intermediate files"
    rm -f "${BAM_PREFIX}.sam" "${BAM_PREFIX}_mapped.sam" "${BAM_PREFIX}_mapped.bam"

    # ---------- Step 4: Filter uniquely mapped reads (optional) ----------
    if [[ "$KEEP_NON_UNIQUE" = true ]]; then
        log_step "Step 4: Keeping all mapped reads (including non-unique)"
        run_command "cp '${BAM_PREFIX}_mapped_sorted.bam' '$FINAL_BAM'" "Copy all mapped reads as final BAM"
        run_command "cp '${BAM_PREFIX}_mapped_sorted.bam.bai' '$FINAL_BAM_INDEX'" "Copy BAM index"
        log_success "Kept all mapped reads (including non-unique)"
    else
        log_step "Step 4: Filtering uniquely mapped reads"
        if [[ -f "../unisRNA.py" ]]; then
            UNSORTED_UNIQUE_BAM="${BAM_PREFIX}_unique_unsorted.bam"
            run_command "python3 '../unisRNA.py' -i '${BAM_PREFIX}_mapped_sorted.bam' -o '$UNSORTED_UNIQUE_BAM'" "Filter unique reads"
            check_file "$UNSORTED_UNIQUE_BAM" "Unsorted unique BAM file"
            run_command "samtools sort -m $SORT_MEMORY --threads $SAMTOOLS_THREADS '$UNSORTED_UNIQUE_BAM' -o '$FINAL_BAM'" "Sort unique BAM file"
            run_command "samtools index -@ $SAMTOOLS_THREADS '$FINAL_BAM'" "Index unique BAM"
            run_command "rm -f '$UNSORTED_UNIQUE_BAM'" "Clean unsorted unique BAM"
            log_success "Filtered to uniquely mapped reads only using unisRNA.py"
        else
            log_warning "unisRNA.py not found, keeping all mapped reads as unique"
            run_command "cp '${BAM_PREFIX}_mapped_sorted.bam' '$FINAL_BAM'" "Copy all mapped reads as final BAM"
            run_command "cp '${BAM_PREFIX}_mapped_sorted.bam.bai' '$FINAL_BAM_INDEX'" "Copy BAM index"
        fi
    fi
	
    # Clean up intermediates
    log_step "Cleaning additional intermediate files"
    rm -f "${BAM_PREFIX}_mapped_sorted.bam" "${BAM_PREFIX}_mapped_sorted.bam.bai"
    rm -f "${INDEX_PREFIX}".*.bt2
    log_success "Upstream analysis completed, BAM file ready for downstream analysis"
fi

# ---------- Step 5: Contig filtering ----------
log_step "Step 5: Filtering contigs based on sRNA mapping"
# Check if filtered contigs already exist and we're not forcing rerun
CONTIG_FILTER_OUTPUT="${CORRELATION_PREFIX}_highly_correlated_contigs.txt"
if [[ "$FORCE_RERUN" = false ]] && [[ -f "$CONTIG_FILTER_OUTPUT" ]]; then
    log_success "Contig filter output file already exists, skipping filtering"
else
    # Build contig_filter arguments
    CONTIG_FILTER_ARGS="--input_bam '$FINAL_BAM' --output '$CORRELATION_PREFIX'"
    CONTIG_FILTER_ARGS="$CONTIG_FILTER_ARGS --min_length $MIN_LENGTH --max_length $MAX_LENGTH"
    CONTIG_FILTER_ARGS="$CONTIG_FILTER_ARGS --threads $CORR_CORES"
    CONTIG_FILTER_ARGS="$CONTIG_FILTER_ARGS --mean_r $MEAN_R --p_value $P_VALUE"
    CONTIG_FILTER_ARGS="$CONTIG_FILTER_ARGS --number_21_nt $NUMBER_21_NT"
    CONTIG_FILTER_ARGS="$CONTIG_FILTER_ARGS --correlation_method $CORRELATION_METHOD"
    CONTIG_FILTER_ARGS="$CONTIG_FILTER_ARGS --mode $VECTOR_MODE"
    
    # Add benchmarks file if provided
    if [[ -n "$BENCHMARK_FILE" ]]; then
        if [[ -f "../$BENCHMARK_FILE" ]]; then
            CONTIG_FILTER_ARGS="$CONTIG_FILTER_ARGS --benchmarks '../$BENCHMARK_FILE'"
        else
            CONTIG_FILTER_ARGS="$CONTIG_FILTER_ARGS --benchmarks '$BENCHMARK_FILE'"
        fi
    fi
    
    run_command "python3 '../contig_filter.py' $CONTIG_FILTER_ARGS" "Contig filtering and correlation analysis"
fi

if [[ ! -f "$CONTIG_FILTER_OUTPUT" ]] || [[ ! -s "$CONTIG_FILTER_OUTPUT" ]]; then
    log_error "CRITICAL: No contigs retained in '$CONTIG_FILTER_OUTPUT'. File exists: $([[ -f "$CONTIG_FILTER_OUTPUT" ]] && echo 'YES' || echo 'NO'), File not empty: $([[ -s "$CONTIG_FILTER_OUTPUT" ]] && echo 'YES' || echo 'NO')"
    exit 1
fi

# ---------- Step 6: Generate feature vectors ----------
log_step "Step 6: Generating ${VECTOR_MODE} feature vectors"
# Check if feature CSV already exists and we're not forcing rerun
if [[ "$FORCE_RERUN" = false ]] && [[ -f "$FEATURE_CSV" ]]; then
    log_success "Feature CSV file already exists, skipping generation"
else
    # Build vector.py arguments
    VECTOR_ARGS="--input '$FINAL_BAM' --output '$FEATURE_CSV'"
    VECTOR_ARGS="$VECTOR_ARGS --min_length $MIN_LENGTH --max_length $MAX_LENGTH"
    VECTOR_ARGS="$VECTOR_ARGS --mode $VECTOR_MODE"
    VECTOR_ARGS="$VECTOR_ARGS --cores $VECTOR_CORES"
    
    # Add filtered contigs file if it exists (from previous step)
    if [[ -f "$CONTIG_FILTER_OUTPUT" ]]; then
        VECTOR_ARGS="$VECTOR_ARGS --filtered_contigs '$CONTIG_FILTER_OUTPUT'"
        log_info "Using filtered contigs file: $CONTIG_FILTER_OUTPUT"
    else
        log_warning "Filtered contigs file not found: $CONTIG_FILTER_OUTPUT, using all contigs from BAM"
    fi
    
    # Add benchmarks file if provided
    if [[ -n "$BENCHMARK_FILE" ]]; then
        if [[ -f "../$BENCHMARK_FILE" ]]; then
            VECTOR_ARGS="$VECTOR_ARGS --benchmarks '../$BENCHMARK_FILE'"
        else
            VECTOR_ARGS="$VECTOR_ARGS --benchmarks '$BENCHMARK_FILE'"
        fi
        log_info "Using benchmark file for feature vector generation"
    fi
    
    run_command "python3 '../vector.py' $VECTOR_ARGS" "Generate ${VECTOR_MODE} feature vectors"
    check_file "$FEATURE_CSV" "Feature CSV file"
fi

# ---------- Step 7: Correlation clustering analysis ----------
log_step "Step 7: ${CORRELATION_DESCRIPTION} clustering analysis"

# ---------- Benchmark file handling ----------
BENCHMARK_ARG=""
if [[ -n "$BENCHMARK_FILE" ]]; then
    if [[ -f "../$BENCHMARK_FILE" ]]; then
        BENCHMARK_ARG="-b ../$BENCHMARK_FILE"
        log_info "Using provided benchmark file from parent directory: ../$BENCHMARK_FILE"
    elif [[ -f "$BENCHMARK_FILE" ]]; then
        BENCHMARK_ARG="-b $BENCHMARK_FILE"
        log_info "Using provided benchmark file from current directory: $BENCHMARK_FILE"
    fi

    # Check file type and log (only if we found a benchmark file)
    if [[ -n "$BENCHMARK_ARG" ]]; then
        if [[ "$BENCHMARK_FILE" == *.csv ]]; then
            log_info "Benchmark file type: CSV format"
        else
            log_info "Benchmark file type: Text file with contig names (one per line)"
        fi
    else
        log_warning "Benchmark file not found: $BENCHMARK_FILE"
        log_info "Proceeding without benchmark file - will use default vsiRNA simulant"
    fi
else
    log_info "No benchmark file specified - will use default vsiRNA simulant"
fi

# Build correlation arguments
CORRELATION_ARGS="-i $FEATURE_CSV $BENCHMARK_ARG --threads $CORR_CORES"
CORRELATION_ARGS="$CORRELATION_ARGS --mean_r $MEAN_R --p_value $P_VALUE"
CORRELATION_ARGS="$CORRELATION_ARGS --min_length $MIN_LENGTH --max_length $MAX_LENGTH"
CORRELATION_ARGS="$CORRELATION_ARGS --correlation_method $CORRELATION_METHOD"
CORRELATION_ARGS="$CORRELATION_ARGS -o '$CORRELATION_PREFIX'"

# Define FILTERED_CONTIGS_FILE (output from corr.py)
FILTERED_CONTIGS_FILE="${CORRELATION_PREFIX}_filtered_contigs.csv"

# Check if clustering results already exist
HEATMAP_FILE="${CORRELATION_PREFIX}_${CORRELATION_METHOD}_clustering_heatmap.svg"
CLUSTER_RESULTS="${CORRELATION_PREFIX}_clustering_order.csv"
FILTERED_CONTIGS_FILE="${CORRELATION_PREFIX}_filtered_contigs.csv"
if [[ "$FORCE_RERUN" = false ]] && [[ -f "$HEATMAP_FILE" ]] && [[ -f "$CLUSTER_RESULTS" ]] && [[ -f "$FILTERED_CONTIGS_FILE" ]]; then
    log_success "Clustering results already exist, skipping clustering analysis"
else
    run_command "python3 '../corr.py' $CORRELATION_ARGS" "${CORRELATION_DESCRIPTION} correlation clustering"
fi

# ---------- Step 8: Enhanced contig feature analysis ----------
log_step "Step 8: Enhanced contig feature analysis"

# Define expected output files
ENHANCED_CONTIGS_FILE="${CORRELATION_PREFIX}_enhanced_contigs.csv"
ENHANCED_SEQUENCES_FILE="${CORRELATION_PREFIX}_enhanced_contigs_sequences.fasta"

# Check if filtered contigs file exists
if [[ -f "$FILTERED_CONTIGS_FILE" ]]; then
    
    # Check if we need to run contig_feature.py
    if [[ "$FORCE_RERUN" = true ]] || [[ ! -f "$ENHANCED_CONTIGS_FILE" ]] || [[ ! -f "$ENHANCED_SEQUENCES_FILE" ]]; then
        log_info "Running enhanced contig feature analysis..."
        log_info "Using filtered contigs file: $FILTERED_CONTIGS_FILE"
        
        # Build contig_feature arguments
        CONTIG_FEATURE_ARGS="--fasta '../$FASTA_FILE' --bam '$FINAL_BAM' --contigs '$FILTERED_CONTIGS_FILE'"
        CONTIG_FEATURE_ARGS="$CONTIG_FEATURE_ARGS --output '$ENHANCED_CONTIGS_FILE' --threads $FEATURE_CORES"
        CONTIG_FEATURE_ARGS="$CONTIG_FEATURE_ARGS --min-length $MIN_LENGTH --max-length $MAX_LENGTH"
        
        # Only add benchmarks-file argument if benchmarks file was provided to the main script
        if [[ -n "$BENCHMARK_FILE" ]]; then
            # Find the actual benchmarks file path
            BENCHMARK_FILE_PATH=""
            if [[ -f "../$BENCHMARK_FILE" ]]; then
                BENCHMARK_FILE_PATH="../$BENCHMARK_FILE"
            elif [[ -f "$BENCHMARK_FILE" ]]; then
                BENCHMARK_FILE_PATH="$BENCHMARK_FILE"
            fi
            
            if [[ -n "$BENCHMARK_FILE_PATH" ]]; then
                CONTIG_FEATURE_ARGS="$CONTIG_FEATURE_ARGS --benchmarks-file '$BENCHMARK_FILE_PATH'"
                log_info "Using benchmarks file: $BENCHMARK_FILE_PATH"
            else
                log_warning "Benchmarks file not found: $BENCHMARK_FILE, skipping benchmarks analysis in contig_feature.py"
            fi
        else
            log_info "No benchmarks file provided to main script, skipping benchmarks analysis in contig_feature.py"
        fi
        
        run_command "python3 '../contig_feature.py' $CONTIG_FEATURE_ARGS" "Enhanced contig feature analysis"
    else
        log_success "Enhanced contig files already exist and --force not set, skipping analysis"
    fi
    
    # Add to keep files (whether newly generated or already existed)
    KEEP_FILES+=("$ENHANCED_CONTIGS_FILE")
    KEEP_FILES+=("${CORRELATION_PREFIX}_enhanced_contigs_summary.txt")
    KEEP_FILES+=("$ENHANCED_SEQUENCES_FILE")
    
    # Also keep benchmarks feature file if it was generated
    if [[ -n "$BENCHMARK_FILE" ]]; then
        BENCHMARKS_FEATURE_FILE="${BENCHMARK_FILE%.*}_feature.csv"
        if [[ -f "$BENCHMARKS_FEATURE_FILE" ]]; then
            KEEP_FILES+=("$BENCHMARKS_FEATURE_FILE")
            log_info "Benchmarks feature file will be kept: $BENCHMARKS_FEATURE_FILE"
        fi
    fi
else
    log_warning "Filtered contigs file not found: $FILTERED_CONTIGS_FILE, skipping enhanced analysis"
fi

# ---------- Clean up and keep only essential files ----------
log_step "Cleaning up - keeping only essential data files"
# List of files to keep (essential data files)
KEEP_FILES=(
    "$FINAL_BAM"
    "$FINAL_BAM_INDEX"
    "$FEATURE_CSV"
    "${CORRELATION_PREFIX}_${CORRELATION_METHOD}_clustering_heatmap.svg"
    "$FILTERED_CONTIGS_FILE"
    "${SAMPLE_PREFIX}_pipeline_parameters.txt"
    "${CORRELATION_PREFIX}_${CORRELATION_METHOD}_distance_matrix.csv"
    "${CORRELATION_PREFIX}_analysis_summary.txt"
    "${CORRELATION_PREFIX}_clustering_order.csv"
    "${CORRELATION_PREFIX}_enhanced_contigs.csv"
    "${CORRELATION_PREFIX}_enhanced_contigs_summary.txt"
    "${BENCHMARK_FILE}_feature.csv"
    "${CORRELATION_PREFIX}_enhanced_contigs_sequences.fasta"
)

# Also keep the intermediate feature file if it was created
if [[ -n "${INTERMEDIATE_FEATURE_CSV+x}" && -f "$INTERMEDIATE_FEATURE_CSV" ]]; then
    KEEP_FILES+=("$INTERMEDIATE_FEATURE_CSV")
fi

# Also keep contig filter output
if [[ -f "$CONTIG_FILTER_OUTPUT" ]]; then
    KEEP_FILES+=("$CONTIG_FILTER_OUTPUT")
    KEEP_FILES+=("${CORRELATION_PREFIX}_detailed_results.csv")
fi

# Also keep all log files (they have timestamps and won't conflict)
KEEP_PATTERNS=(
    "*.log"
)

log_info "Keeping only essential data files:"
for file in "${KEEP_FILES[@]}"; do
    if [[ -f "$file" ]]; then
        log_info "  - $file"
    fi
done

# Remove files that match correlation or sample prefix but aren't in our keep list
find . -maxdepth 1 -type f \( -name "${CORRELATION_PREFIX}_*" -o -name "${SAMPLE_PREFIX}_*" \) | while read -r file; do
    filename=$(basename "$file")
    keep_file=false
    
    # Check if file is in KEEP_FILES
    for keep in "${KEEP_FILES[@]}"; do
        if [[ "$filename" == "$keep" ]]; then
            keep_file=true
            break
        fi
    done
    
    # Check if file matches keep patterns (like log files)
    for pattern in "${KEEP_PATTERNS[@]}"; do
        if [[ "$filename" == $pattern ]]; then
            keep_file=true
            break
        fi
    done
    
    if [[ "$keep_file" == false ]]; then
        log_info "Removing non-essential file: $filename"
        rm -f "$file"
    fi
done

# ---------- Summary ----------
log_step "Pipeline completion report"
log_success "MODULAR sRNA Analysis Pipeline Completed Successfully!"
echo ""
echo "Essential data files retained:"
for file in "${KEEP_FILES[@]}"; do
    if [[ -f "$file" ]]; then
        echo "  - $file"
    fi
done
echo ""
echo "Log files (also retained):"
find . -maxdepth 1 -name "*.log" -type f | while read -r logfile; do
    echo "  - $(basename "$logfile")"
done

# Record resume status in parameters file
RESUME_STATUS="Full analysis run"
if [[ "$SKIP_TO_DOWNSTREAM" = true ]]; then
    RESUME_STATUS="Resumed from existing BAM file"
fi

# Record benchmark method
BENCHMARK_METHOD="Provided benchmark file: $BENCHMARK_FILE"
if [[ -z "$BENCHMARK_FILE" ]]; then
    BENCHMARK_METHOD="Default vsiRNA simulant reference (contigs merged into feature file)"
fi

# Record read filtering method
READ_FILTERING_METHOD="Filtered to uniquely mapped reads only"
if [[ "$KEEP_NON_UNIQUE" = true ]]; then
    READ_FILTERING_METHOD="Kept all mapped reads (including non-unique)"
fi

# Record correlation method
CORRELATION_METHOD_DISPLAY=$CORRELATION_DESCRIPTION

cat > "${SAMPLE_PREFIX}_pipeline_parameters.txt" << EOF
MODULAR sRNA Analysis Pipeline Runtime Parameters
Run time: $(date)
Reference genome (pre-assembled contigs): $FASTA_FILE ($FASTA_BASENAME)
Sequencing data: $READS_FILE ($READS_BASENAME)
Benchmark method: $BENCHMARK_METHOD
Read filtering: $READ_FILTERING_METHOD
Correlation method: $CORRELATION_METHOD_DISPLAY
Output directory: $OUTPUT_DIR
Feature vector mode: $VECTOR_MODE

Performance configuration:
- Total threads: $TOTAL_THREADS
- Bowtie2 threads: $BOWTIE2_THREADS
- Samtools threads: $SAMTOOLS_THREADS
- Sort memory: $SORT_MEMORY
- Vector generation cores: $VECTOR_CORES
- Correlation analysis cores: $CORR_CORES
- Feature analysis cores: $FEATURE_CORES

Analysis parameters:
- sRNA length range: ${MIN_LENGTH}-${MAX_LENGTH}nt
- 21nt threshold: $NUMBER_21_NT
- Mean correlation threshold: $MEAN_R
- P-value threshold: $P_VALUE

Run type: $RESUME_STATUS

Essential output files:
- Final BAM: $FINAL_BAM
- BAM index: $FINAL_BAM_INDEX
- ${VECTOR_MODE} features: $FEATURE_CSV
EOF

# Add intermediate file info if it exists
if [[ -n "${INTERMEDIATE_FEATURE_CSV+x}" && -f "$INTERMEDIATE_FEATURE_CSV" ]]; then
    cat >> "${SAMPLE_PREFIX}_pipeline_parameters.txt" << EOF
- Intermediate features with vsiRNA simulant: $INTERMEDIATE_FEATURE_CSV
EOF
fi

# Add contig filter output info
if [[ -f "$CONTIG_FILTER_OUTPUT" ]]; then
    cat >> "${SAMPLE_PREFIX}_pipeline_parameters.txt" << EOF
- Contig filter output: $CONTIG_FILTER_OUTPUT
EOF
fi

cat >> "${SAMPLE_PREFIX}_pipeline_parameters.txt" << EOF
- Clustering heatmap: ${CORRELATION_PREFIX}_${CORRELATION_METHOD}_clustering_heatmap.svg
- Filtered contigs for correlation: ${CORRELATION_PREFIX}_filtered_contigs.csv
- Enhanced contig features: ${CORRELATION_PREFIX}_enhanced_contigs.csv
- Enhanced contig sequences: ${CORRELATION_PREFIX}_enhanced_contigs_sequences.fasta
- sRNA length range for feature analysis: ${MIN_LENGTH}-${MAX_LENGTH}nt
EOF

log_success "Runtime parameters saved to: ${SAMPLE_PREFIX}_pipeline_parameters.txt"
log_success "Detailed log saved to: $LOG_FILE"