#!/bin/bash
# ==============================================================================
# Script: download_demo.sh
# Description: Downloads demo dataset (HCC1395) for testing the pipeline.
# ==============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

log() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] [DEMO] $*"; }

mkdir -p "${READS_DIR}"

log "Downloading Demo Normal Sample (${NORMAL_ID})..."
# Logic to calculate subdirectory (e.g. 008)
LAST_DIGIT="${NORMAL_ID: -1}"
wget -nc -P "${READS_DIR}" "http://ftp.sra.ebi.ac.uk/vol1/fastq/${NORMAL_ID:0:6}/00${LAST_DIGIT}/${NORMAL_ID}/${NORMAL_ID}_1.fastq.gz"
wget -nc -P "${READS_DIR}" "http://ftp.sra.ebi.ac.uk/vol1/fastq/${NORMAL_ID:0:6}/00${LAST_DIGIT}/${NORMAL_ID}/${NORMAL_ID}_2.fastq.gz"

log "Downloading Demo Tumor Sample (${TUMOR_ID})..."
LAST_DIGIT="${TUMOR_ID: -1}"
wget -nc -P "${READS_DIR}" "http://ftp.sra.ebi.ac.uk/vol1/fastq/${TUMOR_ID:0:6}/00${LAST_DIGIT}/${TUMOR_ID}/${TUMOR_ID}_1.fastq.gz"
wget -nc -P "${READS_DIR}" "http://ftp.sra.ebi.ac.uk/vol1/fastq/${TUMOR_ID:0:6}/00${LAST_DIGIT}/${TUMOR_ID}/${TUMOR_ID}_2.fastq.gz"

# Subsetting for speed testing
log "Subsetting to 5M reads..."
for ID in "$NORMAL_ID" "$TUMOR_ID"; do
    if [ ! -f "${READS_DIR}/${ID}_subset_1.fastq.gz" ]; then
        seqtk sample -s100 "${READS_DIR}/${ID}_1.fastq.gz" 5000000 | gzip > "${READS_DIR}/${ID}_subset_1.fastq.gz"
        seqtk sample -s100 "${READS_DIR}/${ID}_2.fastq.gz" 5000000 | gzip > "${READS_DIR}/${ID}_subset_2.fastq.gz"
    fi
done

log "Demo data ready."
