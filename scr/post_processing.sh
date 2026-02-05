#!/bin/bash

# ==============================================================================
# Script: post_processing.sh
# Description: Annotation and Filtering.
# ==============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] [POST] $*"; }

INPUT_VCF="${RESULTS_DIR}/somatic_filtered.vcf.gz"
FINAL_VCF="${RESULTS_DIR}/somatic_functional.vcf"
FINAL_MAF="${RESULTS_DIR}/somatic_functional.maf"
FOCUS_TYPES="Nonsense_Mutation|Frame_Shift_Del|Frame_Shift_Ins|Splice_Site|Missense_Mutation|In_Frame_Del|In_Frame_Ins|Start_Codon|Nonstop_Mutation"

# ----------------------- [Step 1] Funcotator Prep ---------------------------- #
log "Step 1: Checking Data Sources..."
mkdir -p "${FUNC_DIR}"
FUNC_DATA="${FUNC_DIR}/funcotator_dataSources.v1.8.hg38.20230908s"

if [ ! -d "${FUNC_DATA}" ]; then
    log "Downloading Funcotator data sources..."
    wget -c -P "${FUNC_DIR}" https://storage.googleapis.com/gcp-public-data--broad-references/funcotator/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz
    tar -xzf "${FUNC_DIR}/funcotator_dataSources.v1.8.hg38.20230908s.tar.gz" -C "${FUNC_DIR}"
fi

GTF="${FUNC_DATA}/gencode/hg38/gencode.v43.annotation.REORDERED.gtf"
[ ! -f "${GTF}.idx" ] && gatk IndexFeatureFile -I "${GTF}"

# ----------------------- [Step 2] Annotation --------------------------------- #
log "Step 2: Annotation..."

# VCF Output
if [ ! -f "${RESULTS_DIR}/somatic_annotated.vcf" ]; then
    gatk Funcotator \
        --variant "${INPUT_VCF}" --reference "${REF_FASTA}" --ref-version hg38 \
        --data-sources-path "${FUNC_DATA}" \
        --output "${RESULTS_DIR}/somatic_annotated.vcf" --output-file-format VCF
fi

# MAF Output
if [ ! -f "${RESULTS_DIR}/somatic_full.maf" ]; then
    gatk Funcotator \
        --variant "${INPUT_VCF}" --reference "${REF_FASTA}" --ref-version hg38 \
        --data-sources-path "${FUNC_DATA}" \
        --output "${RESULTS_DIR}/somatic_full.maf" --output-file-format MAF
fi

# ----------------------- [Step 3] Functional Filtering ----------------------- #
log "Step 3: Filtering for high-priority functional variants..."

# Filter VCF (Extracts variants where FUNCOTATION matches the focus types)
grep "^#" "${RESULTS_DIR}/somatic_annotated.vcf" > "${FINAL_VCF}"
grep -Ei "FUNCOTATION=\[[^]]*\|[^]]*\|[^]]*\|[^]]*\|[^]]*\|($FOCUS_TYPES)" "${RESULTS_DIR}/somatic_annotated.vcf" >> "${FINAL_VCF}"

# Filter MAF (Filters Column 9: Variant_Classification)
{
  grep "^#" "${RESULTS_DIR}/somatic_full.maf"
  grep -m 1 "^Hugo_Symbol" "${RESULTS_DIR}/somatic_full.maf"
  awk -F'\t' -v types="$FOCUS_TYPES" 'BEGIN{IGNORECASE=1}
    $0 ~ /^#/ || $0 ~ /^Hugo_Symbol/ {next}
    $9 ~ types {print}
  ' "${RESULTS_DIR}/somatic_full.maf"
} > "${FINAL_MAF}"

# ----------------------- [Step 4] Reporting ---------------------------------- #
log "--------------------------------------------------"
log " [REPORT] Variant counts before filtering (Total):"
grep -v "^#" "${RESULTS_DIR}/somatic_full.maf" | cut -f 9 | sort | uniq -c | sort -nr

log "--------------------------------------------------"
log " [REPORT] Variant counts after filtering (Functional):"
grep -v "^#" "${FINAL_MAF}" | cut -f 9 | sort | uniq -c | sort -nr

log "--------------------------------------------------"
log "[SUCCESS] Post-processing complete."
log "Filtered VCF: ${FINAL_VCF}"
log "Filtered MAF: ${FINAL_MAF}"
