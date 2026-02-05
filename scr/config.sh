#!/bin/bash
# ==============================================================================
# Config: Pipeline Configuration
# Description: Defines paths, sample IDs, and resource URLs.
# ==============================================================================

# --- Project Base Path ---
# Automatically detects the project root (assuming config.sh is in scr/)
# If you run from the root, $(pwd) works.
# If you want to be safe, set this manually or use dynamic resolution:
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

# --- Directory Structure ---
READS_DIR="${PROJECT_DIR}/reads"
RES_DIR="${PROJECT_DIR}/resources"
ALIGNED_DIR="${PROJECT_DIR}/aligned_reads"
RESULTS_DIR="${PROJECT_DIR}/results"
TMP_DIR="${PROJECT_DIR}/tmp"
PLOTS_DIR="${PROJECT_DIR}/plots"

# Sub-directories
RAW_QC_DIR="${READS_DIR}/raw_fastqc"
TRIMMED_DIR="${READS_DIR}/trimmed"
REF_DIR="${RES_DIR}/ref_fasta"
BQSR_RES="${RES_DIR}/BQSR"
MUTECT_RES="${RES_DIR}/mutect2_supporting_files"
FUNC_DIR="${RES_DIR}/funcotator_data_sources"
QC_METRICS_DIR="${ALIGNED_DIR}/post_alignment_qc"

# --- Sample IDs ---
# Users should change these to match their fastq filenames: ${ID}_1.fastq.gz
NORMAL_ID="SRR7890918"
TUMOR_ID="SRR7890919"

# --- Resources & References ---
REF_FASTA="${REF_DIR}/hg38.fa"
GNOMAD_VCF="${MUTECT_RES}/af-only-gnomad.hg38.vcf.gz"
PON_VCF="${MUTECT_RES}/1000g_pon.hg38.vcf.gz"
INTERVALS="${MUTECT_RES}/exome_calling_regions.v1.1.interval_list"

# --- System ---
THREADS=8
