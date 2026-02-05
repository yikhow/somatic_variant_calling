#!/bin/bash
# ==============================================================================
# Script: preprocessing.sh
# Description: Upstream NGS pipeline (QC -> Alignment -> BQSR).
# ==============================================================================

set -euo pipefail
# Load Configuration
source "$(dirname "$0")/config.sh"

ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] [PREP] $*"; }

# Create Directories
mkdir -p "${RAW_QC_DIR}" "${TRIMMED_DIR}" "${ALIGNED_DIR}" \
         "${ALIGNED_DIR}/dedup_bams" "${QC_METRICS_DIR}" \
         "${REF_DIR}" "${BQSR_RES}" "${TMP_DIR}" 

# Define Sample Array for Looping
declare -A SAMPLES=( ["${NORMAL_ID}"]="normal" ["${TUMOR_ID}"]="tumor" )

# ----------------------- [Step 1] Resource Check ----------------------------- #
log "Step 1: Checking Reference and Resources..."

# Download Ref if missing
if [ ! -f "${REF_FASTA}" ]; then
    log "Downloading hg38 reference..."
    wget -nc -P "${REF_DIR}" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip -f "${REF_DIR}/hg38.fa.gz"
fi

# Build Indices
[ ! -f "${REF_FASTA}.fai" ] && samtools faidx "${REF_FASTA}"
[ ! -f "${REF_DIR}/hg38.dict" ] && gatk CreateSequenceDictionary R="${REF_FASTA}" O="${REF_DIR}/hg38.dict"
[ ! -f "${REF_FASTA}.bwt" ] && bwa index "${REF_FASTA}"

# Download BQSR Sites
KNOWN_SITES=(
    "1000G_phase1.snps.high_confidence.hg38.vcf.gz"
    "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    "Homo_sapiens_assembly38.dbsnp138.vcf"
)
for SITE in "${KNOWN_SITES[@]}"; do
    if [ ! -f "${BQSR_RES}/${SITE}" ]; then
        log "Downloading ${SITE}..."
        # Note: Using gsutil requires it to be installed. Alternatively use HTTP mirrors if gsutil is absent.
        gsutil -m cp "gs://gcp-public-data--broad-references/hg38/v0/${SITE}*" "${BQSR_RES}/"
    fi
done

# ----------------------- [Step 2] Raw QC ------------------------------------- #
log "Step 2: Initial FastQC..."
for ID in "${!SAMPLES[@]}"; do
    if [ ! -f "${RAW_QC_DIR}/${ID}_subset_1_fastqc.html" ]; then
        log "Running FastQC for ${ID}..."
        # Using _subset_ for demo flow. Change to original for production.
        fastqc -t "${THREADS}" -o "${RAW_QC_DIR}" "${READS_DIR}/${ID}_subset_1.fastq.gz" "${READS_DIR}/${ID}_subset_2.fastq.gz"
    fi
done

# ----------------------- [Step 3] Trimming ----------------------------------- #
log "Step 3: Adapter Trimming..."
for ID in "${!SAMPLES[@]}"; do
    if [ ! -f "${TRIMMED_DIR}/${ID}_1_trimmed.fastq.gz" ]; then
        log "Trimming ${ID}..."
        fastp -i "${READS_DIR}/${ID}_subset_1.fastq.gz" -I "${READS_DIR}/${ID}_subset_2.fastq.gz" \
              -o "${TRIMMED_DIR}/${ID}_1_trimmed.fastq.gz" -O "${TRIMMED_DIR}/${ID}_2_trimmed.fastq.gz" \
              --html "${TRIMMED_DIR}/${ID}_fastp.html" --json "${TRIMMED_DIR}/${ID}_fastp.json" \
              --thread "${THREADS}"
    fi
done

# ----------------------- [Step 4] Alignment ---------------------------------- #
log "Step 4: BWA Alignment..."
for ID in "${!SAMPLES[@]}"; do
    NAME=${SAMPLES[$ID]}
    if [ ! -f "${ALIGNED_DIR}/${NAME}_sorted.bam" ]; then
        log "Aligning ${NAME}..."
        # SM tag must match what Mutect2 expects (Normal/Tumor) or be handled carefully.
        # Here we use the generic name (normal/tumor) capitalized for SM.
        bwa mem -t "${THREADS}" -M -R "@RG\tID:${ID}\tPL:ILLUMINA\tLB:LIB1\tSM:${NAME^^}" \
            "${REF_FASTA}" "${TRIMMED_DIR}/${ID}_1_trimmed.fastq.gz" "${TRIMMED_DIR}/${ID}_2_trimmed.fastq.gz" \
            | samtools sort -@ 4 -o "${ALIGNED_DIR}/${NAME}_sorted.bam"
        samtools index "${ALIGNED_DIR}/${NAME}_sorted.bam"
    fi
done

# ----------------------- [Step 5] Mark Duplicates ---------------------------- #
log "Step 5: Mark Duplicates..."
for NAME in "${SAMPLES[@]}"; do
    if [ ! -f "${ALIGNED_DIR}/dedup_bams/${NAME}_dedup.bam" ]; then
        gatk MarkDuplicates \
            -I "${ALIGNED_DIR}/${NAME}_sorted.bam" \
            -O "${ALIGNED_DIR}/dedup_bams/${NAME}_dedup.bam" \
            -M "${ALIGNED_DIR}/dedup_bams/${NAME}_metrics.txt" \
            --CREATE_INDEX true
    fi
done

# ----------------------- [Step 6] BQSR --------------------------------------- #
log "Step 6: BQSR..."
mkdir -p "${ALIGNED_DIR}/bqsr_final"
for NAME in "${SAMPLES[@]}"; do
    if [ ! -f "${ALIGNED_DIR}/bqsr_final/${NAME}_final.bam" ]; then
        # Recalibrate
        gatk BaseRecalibrator \
            -I "${ALIGNED_DIR}/dedup_bams/${NAME}_dedup.bam" -R "${REF_FASTA}" \
            --known-sites "${BQSR_RES}/Homo_sapiens_assembly38.dbsnp138.vcf" \
            --known-sites "${BQSR_RES}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
            --known-sites "${BQSR_RES}/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
            -O "${ALIGNED_DIR}/bqsr_final/${NAME}_recal.table"
        
        # Apply
        gatk ApplyBQSR \
            -I "${ALIGNED_DIR}/dedup_bams/${NAME}_dedup.bam" -R "${REF_FASTA}" \
            --bqsr-recal-file "${ALIGNED_DIR}/bqsr_final/${NAME}_recal.table" \
            -O "${ALIGNED_DIR}/bqsr_final/${NAME}_final.bam"
    fi
done

# ----------------------- [Step 7] Final QC ----------------------------------- #
log "Step 7: Final Metrics & MultiQC..."
for NAME in "${SAMPLES[@]}"; do
    [ ! -f "${QC_METRICS_DIR}/${NAME}_alignment_metrics.txt" ] && \
        gatk CollectAlignmentSummaryMetrics -R "${REF_FASTA}" \
        -I "${ALIGNED_DIR}/bqsr_final/${NAME}_final.bam" \
        -O "${QC_METRICS_DIR}/${NAME}_alignment_metrics.txt"
    
    [ ! -f "${QC_METRICS_DIR}/${NAME}_insert_size_metrics.txt" ] && \
        gatk CollectInsertSizeMetrics \
        -I "${ALIGNED_DIR}/bqsr_final/${NAME}_final.bam" \
        -O "${QC_METRICS_DIR}/${NAME}_insert_size_metrics.txt" \
        -H "${QC_METRICS_DIR}/${NAME}_insert_size_histogram.pdf"
done

if [ ! -d "${QC_METRICS_DIR}/multiqc_report" ]; then
    multiqc "${QC_METRICS_DIR}" -o "${QC_METRICS_DIR}/multiqc_report"
fi

log "[SUCCESS] Preprocessing Done."
