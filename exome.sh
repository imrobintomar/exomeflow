#!/bin/bash

################################################################################
# Whole Exome Sequencing Analysis Pipeline
# Author: Robin Tomar, AIIMS New Delhi
# Description: Complete NGS analysis workflow with PARALLEL processing
# Features: Parallel sample processing, checkpointing, resume capability
################################################################################

set -euo pipefail  # Exit on error, undefined variables, and pipe failures
export JAVA_OPTS="-Xmx80g"
export PATH="/media/drprabudh/m1/gatk-4.6.2.0:$PATH"

# ============================================================================
# CONFIGURATION
# ============================================================================

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly TIMESTAMP=$(date +'%Y%m%d_%H%M%S')

# Directory paths
readonly INPUT_DIR="/media/drprabudh/m2/ARM-PRJNA1344179/fastq"
readonly QC_DIR="${INPUT_DIR}/QC"
readonly FASTP_DIR="${INPUT_DIR}/filtered_fastp"
readonly MAP_DIR="${INPUT_DIR}/Mapsam"
readonly VCF_DIR="${INPUT_DIR}/VCF"
readonly LOG_DIR="/media/drprabudh/m2/ARM-PRJNA1344179"
readonly CHECKPOINT_DIR="${LOG_DIR}/.checkpoints"

# Reference files
readonly REF_GENOME="/media/drprabudh/m1/hg38/hg38.fa"
readonly DBSNP="/media/drprabudh/m1/vcf_file/dbsnp.vcf.gz"
readonly INDELS_MILLS="/media/drprabudh/m1/vcf_file/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
readonly INDELS_KNOWN="/media/drprabudh/m1/vcf_file/Homo_sapiens_assembly38.known_indels.vcf.gz"

# Exome capture intervals — pick the BED matching your capture kit:
#   Illumina : /media/drprabudh/m1/Downloads/Illumina_Exome_TargetedRegions_v1.2.hg38.bed
#   Twist Core: /media/drprabudh/m1/Downloads/Twist_Exome_Core_Covered_Targets_hg38.bed
#   Twist Comp: /media/drprabudh/m1/Downloads/Twist_Comprehensive_Exome_Covered_Targets_hg38.bed
readonly EXOME_INTERVALS="/media/drprabudh/m1/Downloads/Illumina_Exome_TargetedRegions_v1.2.hg38.bed"
readonly INTERVAL_PADDING=100

# ANNOVAR paths
readonly ANNOVAR_DIR="/media/drprabudh/m1/annovar"
readonly ANNOVAR_HUMANDB="${ANNOVAR_DIR}/hg38_humandb"

# Log files
readonly LOG_FILE="${LOG_DIR}/analysis_${TIMESTAMP}.log"
readonly ERROR_LOG="${LOG_DIR}/errors_${TIMESTAMP}.log"
readonly PROGRESS_FILE="${LOG_DIR}/progress_${TIMESTAMP}.txt"   # reserved for future use

# Performance settings - ADJUST FOR YOUR SYSTEM
readonly MAX_PARALLEL_JOBS=1          # Number of samples to process in parallel
readonly BWA_THREADS=24               # Threads per BWA job (48 total / 2 jobs = 24 each)
readonly GATK_THREADS=24              # Threads per GATK job
readonly FASTP_THREADS=8              # Threads per fastp job
readonly ANNOVAR_THREADS=24           # Threads per ANNOVAR job

# ANNOVAR annotation databases
readonly ANNOVAR_PROTOCOLS="refGene,dbnsfp47a,clinvar_20240416,gnomad41_exome,gnomad41_genome,avsnp150,cosmic84_coding,exac03"
readonly ANNOVAR_OPERATIONS="g,f,f,f,f,f,f,f"

# ============================================================================
# LOGGING FUNCTIONS
# ============================================================================

log() {
    local level="$1"
    shift
    local message="$@"
    local timestamp=$(date +'%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] [${level}] ${message}" | tee -a "${LOG_FILE}"
}

log_info() {
    log "INFO" "$@"
}

log_warn() {
    log "WARN" "$@"
}

log_error() {
    log "ERROR" "$@"
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: $@" >> "${ERROR_LOG}"
}

log_success() {
    log "SUCCESS" "$@"
}

# ============================================================================
# CHECKPOINT FUNCTIONS
# ============================================================================

init_checkpoints() {
    mkdir -p "$CHECKPOINT_DIR"
    log_info "Checkpoint directory: $CHECKPOINT_DIR"
}

mark_step_complete() {
    local sample="$1"
    local step="$2"
    local checkpoint_file="${CHECKPOINT_DIR}/${sample}.${step}.done"
    touch "$checkpoint_file"
}

is_step_complete() {
    local sample="$1"
    local step="$2"
    local checkpoint_file="${CHECKPOINT_DIR}/${sample}.${step}.done"
    [ -f "$checkpoint_file" ]
}

mark_sample_complete() {
    local sample="$1"
    mark_step_complete "$sample" "COMPLETE"
}

is_sample_complete() {
    local sample="$1"
    is_step_complete "$sample" "COMPLETE"
}

# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

check_dependencies() {
    log_info "Checking dependencies..."
    
    local dependencies=("bwa" "samtools" "gatk" "fastp")
    local missing=()
    
    for cmd in "${dependencies[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            missing+=("$cmd")
        fi
    done
    
    if [ ${#missing[@]} -gt 0 ]; then
        log_error "Missing dependencies: ${missing[*]}"
        exit 1
    fi
    
    log_success "All dependencies found"
}

check_reference_files() {
    log_info "Validating reference files..."
    
    local files=("$REF_GENOME" "$DBSNP" "$INDELS_MILLS" "$INDELS_KNOWN" "$ANNOVAR_HUMANDB")
    
    for file in "${files[@]}"; do
        if [ ! -f "$file" ] && [ ! -d "$file" ]; then
            log_error "Reference file not found: $file"
            exit 1
        fi
    done
    
    log_success "All reference files validated"
}

check_input_files() {
    log_info "Checking input FASTQ files..."
    
    if [ ! -d "$INPUT_DIR" ]; then
        log_error "Input directory not found: $INPUT_DIR"
        exit 1
    fi
    
    local fastq_count=$(find "$INPUT_DIR" -maxdepth 1 -name "*.fastq.gz" | wc -l)
    
    if [ "$fastq_count" -eq 0 ]; then
        log_error "No FASTQ files found in $INPUT_DIR"
        exit 1
    fi
    
    log_info "Found $fastq_count FASTQ files"
}

# ============================================================================
# DIRECTORY SETUP
# ============================================================================

setup_directories() {
    log_info "Setting up output directories..."
    
    local dirs=("$QC_DIR" "$FASTP_DIR" "$MAP_DIR" "$VCF_DIR")
    
    for dir in "${dirs[@]}"; do
        if [ ! -d "$dir" ]; then
            mkdir -p "$dir" || {
                log_error "Failed to create directory: $dir"
                exit 1
            }
        fi
    done
    
    log_success "Directories created/verified"
}

# ============================================================================
# SAMPLE IDENTIFICATION
# ============================================================================

get_sample_names() {
    local samples=()
    
    for f in "$INPUT_DIR"/*.fastq.gz; do
        local basename=$(basename "$f")
        local sample=$(echo "$basename" | sed 's/_[12]\.fastq\.gz$//g')

        if [[ ! " ${samples[@]} " =~ " ${sample} " ]]; then
            samples+=("$sample")
        fi
    done
    
    IFS=$'\n' samples=($(sort <<<"${samples[*]}"))
    unset IFS
    
    printf '%s\n' "${samples[@]}"
}

# ============================================================================
# QUALITY CONTROL
# ============================================================================

run_fastp() {
    local sample="$1"
    
    if is_step_complete "$sample" "fastp"; then
        log_info "fastp already completed for $sample, skipping..."
        return 0
    fi
    
    local r1="${INPUT_DIR}/${sample}_1.fastq.gz"
    local r2="${INPUT_DIR}/${sample}_2.fastq.gz"
    
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        log_error "FASTQ files not found for sample: $sample"
        return 1
    fi
    
    log_info "Running fastp for $sample..."
    
    fastp \
        -i "$r1" \
        -I "$r2" \
        --length_required 50 \
        --qualified_quality_phred 30 \
        -o "${FASTP_DIR}/${sample}_1_filtered.fastq.gz" \
        -O "${FASTP_DIR}/${sample}_2_filtered.fastq.gz" \
        --html "${FASTP_DIR}/${sample}_fastp.html" \
        --json "${FASTP_DIR}/${sample}_fastp.json" \
        -w "$FASTP_THREADS" || {
        log_error "fastp failed for $sample"
        return 1
    }
    
    mark_step_complete "$sample" "fastp"
    log_success "fastp completed for $sample"
    return 0
}

# ============================================================================
# ALIGNMENT
# ============================================================================

run_bwa_mem() {
    local sample="$1"
    
    if is_step_complete "$sample" "bwa"; then
        log_info "BWA MEM already completed for $sample, skipping..."
        return 0
    fi
    
    local r1="${FASTP_DIR}/${sample}_1_filtered.fastq.gz"
    local r2="${FASTP_DIR}/${sample}_2_filtered.fastq.gz"
    local output="${MAP_DIR}/${sample}.bam"
    
    log_info "Running BWA MEM for $sample..."
    
    bwa mem \
        -Y \
        -K 100000000 \
        -t "$BWA_THREADS" \
        -R "@RG\tID:${sample}\tPU:${sample}\tSM:${sample}\tLB:${sample}\tPL:illumina" \
        "$REF_GENOME" \
        "$r1" \
        "$r2" | \
    samtools view \
        -Shb \
        -o "$output" \
        - || {
        log_error "BWA MEM failed for $sample"
        return 1
    }
    
    mark_step_complete "$sample" "bwa"
    log_success "BWA MEM and conversion to BAM completed for $sample"
    return 0
}

# ============================================================================
# BAM PROCESSING
# ============================================================================

sort_bam() {
    local sample="$1"
    
    if is_step_complete "$sample" "sort"; then
        log_info "BAM sorting already completed for $sample, skipping..."
        return 0
    fi
    
    local input="${MAP_DIR}/${sample}.bam"
    local output="${MAP_DIR}/${sample}_sorted.bam"
    
    log_info "Sorting BAM for $sample..."
    
    gatk SortSam \
        -I "$input" \
        -O "$output" \
        -SO coordinate || {
        log_error "SortSam failed for $sample"
        return 1
    }
    
    rm -f "$input"
    mark_step_complete "$sample" "sort"
    log_success "BAM sorting completed for $sample"
    return 0
}

generate_flagstat() {
    local sample="$1"
    
    if is_step_complete "$sample" "flagstat"; then
        log_info "Flagstat already generated for $sample, skipping..."
        return 0
    fi
    
    local bam="${MAP_DIR}/${sample}_sorted.bam"
    local output="${MAP_DIR}/${sample}_flagstat.txt"
    
    log_info "Generating flagstat for $sample..."
    
    samtools flagstat "$bam" > "$output" || {
        log_error "Flagstat generation failed for $sample"
        return 1
    }
    
    mark_step_complete "$sample" "flagstat"
    log_success "Flagstat completed for $sample"
    return 0
}

mark_duplicates() {
    local sample="$1"
    
    if is_step_complete "$sample" "markdup"; then
        log_info "MarkDuplicates already completed for $sample, skipping..."
        return 0
    fi
    
    local input="${MAP_DIR}/${sample}_sorted.bam"
    local output="${MAP_DIR}/${sample}_markdup.bam"
    local metrics="${MAP_DIR}/${sample}_markdup_metrics.txt"
    
    log_info "Marking duplicates for $sample..."
    
    gatk MarkDuplicates \
        -I "$input" \
        -O "$output" \
        -M "$metrics" || {
        log_error "MarkDuplicates failed for $sample"
        return 1
    }
    
    mark_step_complete "$sample" "markdup"
    log_success "MarkDuplicates completed for $sample"
    return 0
}

build_bam_index() {
    local sample="$1"
    
    if is_step_complete "$sample" "index"; then
        log_info "BAM index already created for $sample, skipping..."
        return 0
    fi
    
    local bam="${MAP_DIR}/${sample}_markdup.bam"
    
    log_info "Building BAM index for $sample..."
    
    gatk BuildBamIndex -I "$bam" || {
        log_error "BuildBamIndex failed for $sample"
        return 1
    }
    
    mark_step_complete "$sample" "index"
    log_success "BAM index created for $sample"
    return 0
}

# ============================================================================
# BASE QUALITY SCORE RECALIBRATION
# ============================================================================

run_base_recalibration() {
    local sample="$1"
    
    if is_step_complete "$sample" "bqsr"; then
        log_info "BQSR already completed for $sample, skipping..."
        return 0
    fi
    
    local bam="${MAP_DIR}/${sample}_markdup.bam"
    local recal_table="${MAP_DIR}/${sample}_recal.table"
    
    log_info "Running BaseRecalibrator for $sample..."
    
    gatk BaseRecalibrator \
        -I "$bam" \
        --known-sites "$DBSNP" \
        --known-sites "$INDELS_KNOWN" \
        --known-sites "$INDELS_MILLS" \
        -O "$recal_table" \
        -R "$REF_GENOME" || {
        log_error "BaseRecalibrator failed for $sample"
        return 1
    }
    
    log_info "Applying BQSR for $sample..."
    
    local input_bam="${MAP_DIR}/${sample}_markdup.bam"
    local output_bam="${MAP_DIR}/${sample}_recalibrated.bam"
    
    gatk ApplyBQSR \
        --bqsr-recal-file "$recal_table" \
        -I "$input_bam" \
        -O "$output_bam" \
        -R "$REF_GENOME" || {
        log_error "ApplyBQSR failed for $sample"
        return 1
    }
    
    # Index the recalibrated BAM
    samtools index "$output_bam" || {
        log_error "Indexing recalibrated BAM failed for $sample"
        return 1
    }

    # Remove intermediate BAMs to free disk space
    # KEPT: ${sample}_recalibrated.bam + .bai  — use this in IGV for variant validation
    rm -f "${MAP_DIR}/${sample}_sorted.bam" \
          "${MAP_DIR}/${sample}_markdup.bam" \
          "${MAP_DIR}/${sample}_markdup.bai" \
          "$recal_table"
    log_info "[$sample] Kept for IGV: ${MAP_DIR}/${sample}_recalibrated.bam"

    mark_step_complete "$sample" "bqsr"
    log_success "BQSR completed for $sample"
    return 0
}

# ============================================================================
# VARIANT CALLING
# ============================================================================

run_haplotype_caller() {
    local sample="$1"

    if is_step_complete "$sample" "haplotype"; then
        log_info "HaplotypeCaller already completed for $sample, skipping..."
        return 0
    fi

    local bam="${MAP_DIR}/${sample}_recalibrated.bam"
    local vcf="${VCF_DIR}/${sample}.vcf"

    log_info "Running HaplotypeCaller for $sample..."

    # Build interval args if BED file exists
    local interval_args=""
    if [ -f "$EXOME_INTERVALS" ]; then
        interval_args="-L $EXOME_INTERVALS --interval-padding $INTERVAL_PADDING"
        log_info "Using exome intervals: $EXOME_INTERVALS (padding: ${INTERVAL_PADDING}bp)"
    else
        log_warn "No exome intervals BED file found at $EXOME_INTERVALS — calling whole genome (slower, more false positives)"
    fi

    gatk HaplotypeCaller \
        -R "$REF_GENOME" \
        -I "$bam" \
        -O "$vcf" \
        --dbsnp "$DBSNP" \
        $interval_args \
        --native-pair-hmm-threads "$GATK_THREADS" || {
        log_error "HaplotypeCaller failed for $sample"
        return 1
    }

    mark_step_complete "$sample" "haplotype"
    log_success "HaplotypeCaller completed for $sample"
    return 0
}

# ============================================================================
# VARIANT FILTERING (GATK Hard Filters)
# ============================================================================

run_variant_filtration() {
    local sample="$1"

    if is_step_complete "$sample" "filter"; then
        log_info "Variant filtering already completed for $sample, skipping..."
        return 0
    fi

    local vcf="${VCF_DIR}/${sample}.vcf"
    local snp_raw="${VCF_DIR}/${sample}_snp_raw.vcf"
    local indel_raw="${VCF_DIR}/${sample}_indel_raw.vcf"
    local snp_filtered="${VCF_DIR}/${sample}_snp_filtered.vcf"
    local indel_filtered="${VCF_DIR}/${sample}_indel_filtered.vcf"
    local pass_vcf="${VCF_DIR}/${sample}_PASS.vcf"

    # --- Step 1: Separate SNPs and INDELs ---
    log_info "[$sample] Separating SNPs and INDELs..."

    gatk SelectVariants -R "$REF_GENOME" -V "$vcf" \
        --select-type-to-include SNP \
        -O "$snp_raw" || { log_error "SelectVariants (SNP) failed for $sample"; return 1; }

    gatk SelectVariants -R "$REF_GENOME" -V "$vcf" \
        --select-type-to-include INDEL \
        -O "$indel_raw" || { log_error "SelectVariants (INDEL) failed for $sample"; return 1; }

    # --- Step 2: Filter SNPs (GATK best-practice hard filters) ---
    log_info "[$sample] Filtering SNPs..."
    gatk VariantFiltration \
        -R "$REF_GENOME" \
        -V "$snp_raw" \
        -O "$snp_filtered" \
        --filter-expression "QD < 2.0"          --filter-name "SNP_LowQD" \
        --filter-expression "FS > 60.0"          --filter-name "SNP_StrandBias" \
        --filter-expression "SOR > 3.0"          --filter-name "SNP_StrandOddsRatio" \
        --filter-expression "MQ < 40.0"          --filter-name "SNP_LowMQ" \
        --filter-expression "MQRankSum < -12.5"  --filter-name "SNP_MQRankSum" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "SNP_ReadPosRankSum" \
        --filter-expression "DP < 10"            --filter-name "LowDepth" \
        --genotype-filter-expression "GQ < 20"   --genotype-filter-name "LowGQ" || {
        log_error "VariantFiltration (SNP) failed for $sample"; return 1; }

    # --- Step 3: Filter INDELs (different thresholds — no MQ filters) ---
    log_info "[$sample] Filtering INDELs..."
    gatk VariantFiltration \
        -R "$REF_GENOME" \
        -V "$indel_raw" \
        -O "$indel_filtered" \
        --filter-expression "QD < 2.0"             --filter-name "INDEL_LowQD" \
        --filter-expression "FS > 200.0"           --filter-name "INDEL_StrandBias" \
        --filter-expression "SOR > 10.0"           --filter-name "INDEL_StrandOddsRatio" \
        --filter-expression "ReadPosRankSum < -20.0" --filter-name "INDEL_ReadPosRankSum" \
        --filter-expression "DP < 10"              --filter-name "LowDepth" \
        --genotype-filter-expression "GQ < 20"     --genotype-filter-name "LowGQ" || {
        log_error "VariantFiltration (INDEL) failed for $sample"; return 1; }

    # --- Step 4: Merge filtered SNPs + INDELs ---
    log_info "[$sample] Merging filtered SNPs and INDELs..."
    gatk MergeVcfs \
        -I "$snp_filtered" \
        -I "$indel_filtered" \
        -O "${VCF_DIR}/${sample}_merged_filtered.vcf" || {
        log_error "MergeVcfs failed for $sample"; return 1; }

    # --- Step 5: Extract PASS-only variants ---
    log_info "[$sample] Extracting PASS variants..."
    gatk SelectVariants \
        -R "$REF_GENOME" \
        -V "${VCF_DIR}/${sample}_merged_filtered.vcf" \
        -O "$pass_vcf" \
        --exclude-filtered \
        --exclude-non-variants || {
        log_error "SelectVariants (PASS) failed for $sample"; return 1; }

    # --- Step 6: Summary ---
    # grep -vc exits 1 on 0 matches (but still prints "0"), so use || true not || echo 0
    local total; total=$(grep -vc "^#" "$vcf" 2>/dev/null || true); total=${total:-0}
    local snp_count; snp_count=$(grep -vc "^#" "$snp_raw" 2>/dev/null || true); snp_count=${snp_count:-0}
    local indel_count; indel_count=$(grep -vc "^#" "$indel_raw" 2>/dev/null || true); indel_count=${indel_count:-0}
    local passed; passed=$(grep -vc "^#" "$pass_vcf" 2>/dev/null || true); passed=${passed:-0}
    log_info "[$sample] Total: $total | SNPs: $snp_count | INDELs: $indel_count | PASS: $passed | Filtered: $(( total - passed ))"

    # Remove intermediate VCFs and their index files (keep raw, PASS, and annotated only)
    # KEPT: ${sample}.vcf            — raw HaplotypeCaller output (all variants)
    # KEPT: ${sample}_PASS.vcf       — PASS-only variants (fed to ANNOVAR)
    rm -f "$snp_raw" "${snp_raw}.idx" \
          "$indel_raw" "${indel_raw}.idx" \
          "$snp_filtered" "${snp_filtered}.idx" \
          "$indel_filtered" "${indel_filtered}.idx" \
          "${VCF_DIR}/${sample}_merged_filtered.vcf" \
          "${VCF_DIR}/${sample}_merged_filtered.vcf.idx"

    mark_step_complete "$sample" "filter"
    log_success "Variant filtering completed for $sample"
    return 0
}

# ============================================================================
# VARIANT ANNOTATION
# ============================================================================

run_annovar_annotation() {
    local sample="$1"

    if is_step_complete "$sample" "annovar"; then
        log_info "ANNOVAR annotation already completed for $sample, skipping..."
        return 0
    fi

    local vcf="${VCF_DIR}/${sample}_PASS.vcf"
    local output_prefix="${VCF_DIR}/${sample}.annovar"
    
    log_info "Running ANNOVAR annotation for $sample..."
    
    "${ANNOVAR_DIR}/table_annovar.pl" \
        "$vcf" \
        "$ANNOVAR_HUMANDB" \
        --buildver hg38 \
        --out "$output_prefix" \
        --remove \
        --protocol "$ANNOVAR_PROTOCOLS" \
        --operation "$ANNOVAR_OPERATIONS" \
        -nastring . \
        --polish \
        --otherinfo \
        --vcfinput \
        --thread "$ANNOVAR_THREADS" || {
        log_error "ANNOVAR annotation failed for $sample"
        return 1
    }
    
    # Remove ANNOVAR intermediate file (--remove doesn't always clean this up)
    rm -f "${output_prefix}.avinput"

    mark_step_complete "$sample" "annovar"
    log_success "ANNOVAR annotation completed for $sample"
    return 0
}

# ============================================================================
# MAIN ANALYSIS PIPELINE
# ============================================================================

process_sample() {
    local sample="$1"
    local sample_log="${LOG_DIR}/${sample}_${TIMESTAMP}.log"
    
    {
        log_info "=========================================="
        log_info "Processing sample: $sample"
        log_info "=========================================="
        
        run_fastp "$sample" || { log_error "Failed at fastp for $sample"; return 1; }
        run_bwa_mem "$sample" || { log_error "Failed at BWA MEM for $sample"; return 1; }
        sort_bam "$sample" || { log_error "Failed at sort for $sample"; return 1; }
        generate_flagstat "$sample" || { log_error "Failed at flagstat for $sample"; return 1; }
        mark_duplicates "$sample" || { log_error "Failed at markdup for $sample"; return 1; }
        build_bam_index "$sample" || { log_error "Failed at index for $sample"; return 1; }
        run_base_recalibration "$sample" || { log_error "Failed at BQSR for $sample"; return 1; }
        run_haplotype_caller "$sample" || { log_error "Failed at HaplotypeCaller for $sample"; return 1; }
        run_variant_filtration "$sample" || { log_error "Failed at filtering for $sample"; return 1; }
        run_annovar_annotation "$sample" || { log_error "Failed at ANNOVAR for $sample"; return 1; }
        
        mark_sample_complete "$sample"
        log_success "Sample $sample completed successfully"
        log_info "========== OUTPUT FILES: $sample =========="
        log_info "  BAM (IGV):        ${MAP_DIR}/${sample}_recalibrated.bam"
        log_info "  BAM index:        ${MAP_DIR}/${sample}_recalibrated.bam.bai"
        log_info "  Raw VCF:          ${VCF_DIR}/${sample}.vcf"
        log_info "  PASS VCF:         ${VCF_DIR}/${sample}_PASS.vcf"
        log_info "  Annotated VCF:    ${VCF_DIR}/${sample}.annovar.hg38_multianno.vcf"
        log_info "  Annotated TXT:    ${VCF_DIR}/${sample}.annovar.hg38_multianno.txt"
        log_info "==========================================="
        return 0

    } >> "$sample_log" 2>&1
}

main() {
    log_info "=========================================="
    log_info "Whole Exome  Analysis Pipeline"
    log_info "PARALLEL MODE - Processing up to $MAX_PARALLEL_JOBS samples simultaneously"
    log_info "=========================================="
    log_info "Started at: $(date +'%Y-%m-%d %H:%M:%S')"
    log_info "Log file: $LOG_FILE"
    
    # Pre-flight checks
    check_dependencies
    check_reference_files
    check_input_files
    setup_directories
    init_checkpoints
    
    # Get sample names
    log_info "Identifying sample names..."
    local samples=()
    mapfile -t samples < <(get_sample_names)
    
    local total_samples=${#samples[@]}
    
    if [ "$total_samples" -eq 0 ]; then
        log_error "No samples found to process"
        exit 1
    fi
    
    log_info "Found $total_samples sample(s) to process: ${samples[*]}"
    log_info "Will process $MAX_PARALLEL_JOBS sample(s) in parallel"
    
    # Process samples in parallel batches
    local processed=0
    local failed=0
    local active_pids=()
    local active_samples=()
    
    for sample in "${samples[@]}"; do
        # Check if already completed
        if is_sample_complete "$sample"; then
            log_info "Sample $sample already completed, skipping"
            processed=$(( processed + 1 ))
            continue
        fi
        
        # Wait if we've reached max parallel jobs
        if [ ${#active_pids[@]} -ge $MAX_PARALLEL_JOBS ]; then
            log_info "Waiting for a job to complete (${#active_pids[@]} running)..."
            
            # Wait for any job to finish
            wait -n 2>/dev/null || true
            
            # Clean up completed jobs from tracking arrays
            local new_pids=()
            local new_samples=()
            for i in "${!active_pids[@]}"; do
                if kill -0 "${active_pids[$i]}" 2>/dev/null; then
                    new_pids+=("${active_pids[$i]}")
                    new_samples+=("${active_samples[$i]}")
                else
                    # Job finished
                    if is_sample_complete "${active_samples[$i]}"; then
                        log_success "Sample ${active_samples[$i]} completed in background"
                        processed=$(( processed + 1 ))
                    else
                        log_error "Sample ${active_samples[$i]} failed"
                        failed=$(( failed + 1 ))
                    fi
                fi
            done
            active_pids=("${new_pids[@]+"${new_pids[@]}"}")
            active_samples=("${new_samples[@]+"${new_samples[@]}"}")
        fi
        
        # Start new job
        log_info "Starting processing of sample: $sample (${#active_pids[@]}/$MAX_PARALLEL_JOBS active)"
        process_sample "$sample" &
        active_pids+=($!)
        active_samples+=("$sample")
    done
    
    # Wait for remaining jobs
    log_info "Waiting for remaining ${#active_pids[@]} job(s) to complete..."
    for pid in "${active_pids[@]}"; do
        wait "$pid" 2>/dev/null || true
    done
    
    # Check final status
    for sample in "${active_samples[@]}"; do
        if is_sample_complete "$sample"; then
            processed=$(( processed + 1 ))
        else
            failed=$(( failed + 1 ))
        fi
    done
    
    # Summary
    log_info "=========================================="
    log_info "Pipeline Summary"
    log_info "=========================================="
    log_info "Total samples: $total_samples"
    log_info "Successfully processed: $processed"
    log_info "Failed: $failed"
    log_info "Completed at: $(date +'%Y-%m-%d %H:%M:%S')"
    log_info "=========================================="
    
    if [ "$failed" -gt 0 ]; then
        log_error "Pipeline completed with $failed error(s)"
        log_info "Check individual sample logs in: $LOG_DIR"
        exit 1
    else
        log_success "Pipeline completed successfully!"
        log_info "Annotated VCF files are in: $VCF_DIR"
        exit 0
    fi
}

# ============================================================================
# EXECUTION
# ============================================================================

trap 'log_error "Pipeline interrupted"; exit 1' INT TERM

main "$@"
