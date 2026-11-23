#!/usr/bin/env bash

# step2_loop_patients.sh

# Run step2_filter_forcecalled_snvs.R for all patients with a completed Step 1 force-called VCF.

# Usage:

#   bash scripts/02_run_all_patients.sh \
#       FORCECALLED_DIR \
#       SAMPLE_METADATA_TSV \
#       OUTPUT_DIR

# Example:

#   bash scripts/02_run_all_patients.sh \
#       results/01_forcecalled \
#       config/sample_metadata.tsv \
#       results/02_filtered_snvs

# Expected Step 1 structure:

# results/01_forcecalled/
# ├── PAT001/
# │   └── PAT001_mpileup.vcf.gz
# ├── PAT002/
# │   └── PAT002_mpileup.vcf.gz
# └── ...

set -Eeuo pipefail

# 1. Command-line arguments

if [[ $# -ne 3 ]]; then
    echo "Usage:"
    echo "  bash $0 FORCECALLED_DIR SAMPLE_METADATA_TSV OUTPUT_DIR"
    exit 1
fi

FORCECALLED_DIR="$1"
SAMPLE_METADATA="$2"
OUTPUT_DIR="$3"

# 2. Locate the Step 2 R script

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

STEP2_SCRIPT="${SCRIPT_DIR}/02_filter_forcecalled_snvs.R"

# 3. Validate inputs

if [[ ! -d "$FORCECALLED_DIR" ]]; then
    echo "ERROR: Force-called directory not found:"
    echo "  $FORCECALLED_DIR"
    exit 1
fi

if [[ ! -f "$SAMPLE_METADATA" ]]; then
    echo "ERROR: Sample metadata file not found:"
    echo "  $SAMPLE_METADATA"
    exit 1
fi

if [[ ! -f "$STEP2_SCRIPT" ]]; then
    echo "ERROR: Step 2 R script not found:"
    echo "  $STEP2_SCRIPT"
    exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript was not found in PATH."
    exit 1
fi


mkdir -p "$OUTPUT_DIR"

LOG_DIR="${OUTPUT_DIR}/logs"
mkdir -p "$LOG_DIR"

MASTER_LOG="${LOG_DIR}/step2_$(date +%Y%m%d_%H%M%S).log"


# 4. Find all completed Step 1 VCFs

mapfile -t VCF_FILES < <(
    find "$FORCECALLED_DIR" \
        -type f \
        -name "*_mpileup.vcf.gz" \
        | sort
)

if (( ${#VCF_FILES[@]} == 0 )); then
    echo "ERROR: No *_mpileup.vcf.gz files found in:"
    echo "  $FORCECALLED_DIR"
    exit 1
fi

N_PATIENTS=${#VCF_FILES[@]}

echo "" | tee -a "$MASTER_LOG"
echo "Step 2 cohort run"                                           | tee -a "$MASTER_LOG"
echo "Started: $(date)"                                          | tee -a "$MASTER_LOG"
echo "Patients found: $N_PATIENTS"                               | tee -a "$MASTER_LOG"
echo "" | tee -a "$MASTER_LOG"

# 5. Process each patient

n_success=0
n_failed=0
failed_patients=()

i=0

for VCF_PATH in "${VCF_FILES[@]}"; do

    i=$((i + 1))

    VCF_NAME="$(basename "$VCF_PATH")"

    PATIENT="${VCF_NAME%_mpileup.vcf.gz}"

    PATIENT_OUT="${OUTPUT_DIR}/${PATIENT}"
    PATIENT_LOG="${LOG_DIR}/${PATIENT}.log"

    mkdir -p "$PATIENT_OUT"

    echo "" | tee -a "$MASTER_LOG"
    echo "" | tee -a "$MASTER_LOG"
    echo "[$i/$N_PATIENTS] Patient: $PATIENT"                       | tee -a "$MASTER_LOG"
    echo "VCF: $VCF_PATH"                                           | tee -a "$MASTER_LOG"
    echo "Started: $(date)"                                         | tee -a "$MASTER_LOG"


    if Rscript "$STEP2_SCRIPT" \
        "$PATIENT" \
        "$VCF_PATH" \
        "$SAMPLE_METADATA" \
        "$PATIENT_OUT" \
        > "$PATIENT_LOG" 2>&1
    then

        n_success=$((n_success + 1))

        echo "SUCCESS: $PATIENT" | tee -a "$MASTER_LOG"

        tail -10 "$PATIENT_LOG" | tee -a "$MASTER_LOG"

    else

        n_failed=$((n_failed + 1))

        failed_patients+=("$PATIENT")

        echo "FAILED: $PATIENT" | tee -a "$MASTER_LOG"
        echo "Log: $PATIENT_LOG" | tee -a "$MASTER_LOG"

        echo "Last 20 lines:" | tee -a "$MASTER_LOG"

        tail -20 "$PATIENT_LOG" \
            | sed 's/^/    /' \
            | tee -a "$MASTER_LOG"

    fi

done

# 6. Final cohort summary

echo "" | tee -a "$MASTER_LOG"
echo "" | tee -a "$MASTER_LOG"
echo "Step 2 cohort run complete"                                  | tee -a "$MASTER_LOG"
echo "Finished: $(date)"                                          | tee -a "$MASTER_LOG"
echo "" | tee -a "$MASTER_LOG"

echo "Total patients: $N_PATIENTS" | tee -a "$MASTER_LOG"
echo "Succeeded:      $n_success"  | tee -a "$MASTER_LOG"
echo "Failed:         $n_failed"   | tee -a "$MASTER_LOG"


if (( n_failed > 0 )); then

    echo "Failed patients:" | tee -a "$MASTER_LOG"

    printf '  %s\n' "${failed_patients[@]}" \
        | tee -a "$MASTER_LOG"

fi

echo "" | tee -a "$MASTER_LOG"
echo "Master log:" | tee -a "$MASTER_LOG"
echo "  $MASTER_LOG" | tee -a "$MASTER_LOG"

echo "Per-patient logs:" | tee -a "$MASTER_LOG"
echo "  ${LOG_DIR}/<PATIENT>.log" | tee -a "$MASTER_LOG"