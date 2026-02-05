#!/bin/bash
# ==============================================================================
# Script: variant_calling.sh
# Description: Mutect2 Variant Calling.
# ==============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

ts() { date "+%Y-%m-%d %H:%M:%S"; }
log() { echo "[$(ts)] [CALL] $*"; }

mkdir -p "${RESULTS_DIR}" "${MUTECT_RES}" "${TMP_DIR}"

# ----------------------- [Step 0] Prepare Resources -------------------------- #
log "Step 0: Checking Mutect2 Resources..."
RESOURCES=(
    "af-only-gnomad.hg38.vcf.gz|https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
    "af-only-gnomad.hg38.vcf.gz.tbi|https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi"
    "1000g_pon.hg38.vcf.gz|https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
    "1000g_pon.hg38.vcf.gz.tbi|https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi"
    "exome_calling_regions.v1.1.interval_list|https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list"
)

for ENTRY in "${RESOURCES[@]}"; do
    FILE="${ENTRY%%|*}"
    URL="${ENTRY#*|}"
    if [ ! -f "${MUTECT_RES}/${FILE}" ]; then
        log "Downloading ${FILE}..."
        wget -nc -P "${MUTECT_RES}" "${URL}"
    fi
done

# ----------------------- [Step 1] Mutect2 ------------------------------------ #
log "Step 1: Running Mutect2..."

TUMOR_BAM="${ALIGNED_DIR}/bqsr_final/tumor_final.bam"
NORMAL_BAM="${ALIGNED_DIR}/bqsr_final/normal_final.bam"

if [ ! -f "${RESULTS_DIR}/somatic_raw.vcf.gz" ]; then
    gatk Mutect2 \
        -R "${REF_FASTA}" \
        -I "${TUMOR_BAM}" \
        -I "${NORMAL_BAM}" \
        -normal "NORMAL" \
        --germline-resource "${GNOMAD_VCF}" \
        --panel-of-normals "${PON_VCF}" \
        -L "${INTERVALS}" \
        --f1r2-tar-gz "${RESULTS_DIR}/f1r2.tar.gz" \
        -O "${RESULTS_DIR}/somatic_raw.vcf.gz"
fi

# ----------------------- [Step 2] Contamination ------------------------------ #
log "Step 2: Contamination Estimation..."
for SAMPLE in tumor normal; do
    if [ ! -f "${RESULTS_DIR}/${SAMPLE}_pileups.table" ]; then
        BAM_FILE="${ALIGNED_DIR}/bqsr_final/${SAMPLE}_final.bam"
        gatk GetPileupSummaries \
            --java-options '-Xmx8G' --tmp-dir "${TMP_DIR}" \
            -I "${BAM_FILE}" \
            -V "${GNOMAD_VCF}" \
            -L "${INTERVALS}" \
            -O "${RESULTS_DIR}/${SAMPLE}_pileups.table"
    fi
done

if [ ! -f "${RESULTS_DIR}/contamination.table" ]; then
    gatk CalculateContamination \
        -I "${RESULTS_DIR}/tumor_pileups.table" \
        -matched "${RESULTS_DIR}/normal_pileups.table" \
        -O "${RESULTS_DIR}/contamination.table"
fi

# ----------------------- [Step 3] Artifacts & Filtering ---------------------- #
log "Step 3: Filtering..."
if [ ! -f "${RESULTS_DIR}/read-orientation-model.tar.gz" ]; then
    gatk LearnReadOrientationModel \
        -I "${RESULTS_DIR}/f1r2.tar.gz" \
        -O "${RESULTS_DIR}/read-orientation-model.tar.gz"
fi

if [ ! -f "${RESULTS_DIR}/somatic_filtered.vcf.gz" ]; then
    gatk FilterMutectCalls \
        -R "${REF_FASTA}" \
        -V "${RESULTS_DIR}/somatic_raw.vcf.gz" \
        --contamination-table "${RESULTS_DIR}/contamination.table" \
        --ob-priors "${RESULTS_DIR}/read-orientation-model.tar.gz" \
        -O "${RESULTS_DIR}/somatic_filtered.vcf.gz"
fi

log "[SUCCESS] Variant Calling Done."
