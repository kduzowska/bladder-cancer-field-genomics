#!/bin/bash
# Loop over all patients with WES data

BAM_BASE=/path/to/BAM/files
SCRIPTS_DIR=/path/to/scripts
LOG_DIR=$(dirname $SCRIPTS_DIR)/logs
mkdir -p $LOG_DIR

# Master log
MASTER_LOG=$LOG_DIR/master_$(date +%Y%m%d_%H%M%S).log

echo "" | tee -a $MASTER_LOG
echo "Master run started: $(date)" | tee -a $MASTER_LOG
echo "" | tee -a $MASTER_LOG

# Get patient list
# Extract unique patient IDs from BAM folder names
# Pattern: XG-4050-{PATIENT_ID}-{TISSUE}
PATIENTS=$(ls $BAM_BASE | grep "^XG-4050-" | sed 's/^XG-4050-//' | sed 's/-[^-]*$//' | sort -u)

n_patients=$(echo "$PATIENTS" | wc -l)
echo "Found $n_patients unique patients" | tee -a $MASTER_LOG
echo "$PATIENTS" | tee -a $MASTER_LOG
echo "" | tee -a $MASTER_LOG

# Process each patient
i=0
n_success=0
n_failed=0
failed_patients=""

for PATIENT in $PATIENTS; do
    i=$((i+1))
    echo "" | tee -a $MASTER_LOG
    echo "[$i/$n_patients] Patient: $PATIENT — $(date)" | tee -a $MASTER_LOG
    echo "" | tee -a $MASTER_LOG
    
    # Per-patient log
    PATIENT_LOG=$LOG_DIR/${PATIENT}.log
    
    # Run, capturing log
    if bash $SCRIPTS_DIR/step1_snv_calling.sh $PATIENT > $PATIENT_LOG 2>&1; then
        n_success=$((n_success+1))
        echo "  ✓ SUCCESS for $PATIENT" | tee -a $MASTER_LOG
        # Append per-patient summary to master log
        tail -5 $PATIENT_LOG | tee -a $MASTER_LOG
    else
        n_failed=$((n_failed+1))
        failed_patients="$failed_patients $PATIENT"
        echo "  ✗ FAILED for $PATIENT (see $PATIENT_LOG)" | tee -a $MASTER_LOG
        # Show last 10 lines of failed log
        echo "  Last 10 lines of log:" | tee -a $MASTER_LOG
        tail -10 $PATIENT_LOG | sed 's/^/    /' | tee -a $MASTER_LOG
    fi
    echo "" | tee -a $MASTER_LOG
done

# Final summary
echo "" | tee -a $MASTER_LOG
echo "MASTER RUN COMPLETE: $(date)" | tee -a $MASTER_LOG
echo "" | tee -a $MASTER_LOG
echo "Total patients:   $n_patients" | tee -a $MASTER_LOG
echo "Succeeded:        $n_success" | tee -a $MASTER_LOG
echo "Failed:           $n_failed" | tee -a $MASTER_LOG
if [ $n_failed -gt 0 ]; then
    echo "Failed patients: $failed_patients" | tee -a $MASTER_LOG
fi
echo "" | tee -a $MASTER_LOG
echo "Master log: $MASTER_LOG" | tee -a $MASTER_LOG
echo "Per-patient logs: $LOG_DIR/<patient>.log" | tee -a $MASTER_LOG