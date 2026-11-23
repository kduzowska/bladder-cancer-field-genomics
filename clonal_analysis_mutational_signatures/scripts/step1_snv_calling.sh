#!/usr/bin/env bash

# Process one patient:
#   1. Build the union of SNVs identified by Octopus
#   2. Measure the same REF/ALT alleles across all patient BAMs

# Usage:
#   bash step1_process_patients.sh PATIENT_ID

# The vcfs after Octopus calling and annotation have the following names: "PATIENT_*vsBLOOD_somatic.vcf.hg38_multianno.vcf"

# Output:
#   <PATIENT>_union_snvs.vcf.gz - every unique SNV found by Octopus in any tissue vs blood comparison for this patient
#   <PATIENT>_target_alleles.tsv.gz - the variants but in a format for bcftools further (for example: chr1 2742947 A,G)
#   <PATIENT>_mpileup.vcf.gz - for further analyses. Contains each variant with AD, DP

set -Eeuo pipefail # if something suspicious happens, stop instead of quietly continuing

PATIENT="${1:?Usage: $0 PATIENT_ID}"

# Configuration

REF="/path/to/ref/file"

BAM_BASE="/path/to/bam/files"

VCF_BASE="/path/to/vcfs/after/octopus/Output_octopus_ann"

OUT_BASE="/path/to/output/directory"

OUT_DIR="${OUT_BASE}/${PATIENT}"
mkdir -p "$OUT_DIR"

STATUS_FILE="${OUT_DIR}/STATUS.txt"

UNION_VCF="${OUT_DIR}/${PATIENT}_union_snvs.vcf.gz"
ALLELES_TSV="${OUT_DIR}/${PATIENT}_target_alleles.tsv.gz"
OUT_VCF="${OUT_DIR}/${PATIENT}_mpileup.vcf.gz"
ERR_FILE="${OUT_DIR}/${PATIENT}_mpileup.err"

TMP_DIR=""
cleanup() {
    if [[ -n "${TMP_DIR}" && -d "${TMP_DIR}" ]]; then
        rm -rf "${TMP_DIR}"
    fi
}
trap cleanup EXIT

echo "STARTED $(date)" > "$STATUS_FILE"

echo "Patient: ${PATIENT}"
echo "Started: $(date)"
echo "Output:  ${OUT_DIR}"

# Check required files/programs

for cmd in bcftools bgzip tabix; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $cmd"
        echo "FAILED missing_program_${cmd} $(date)" >> "$STATUS_FILE"
        exit 1
    fi
done

if [[ ! -f "$REF" ]]; then
    echo "ERROR: Reference FASTA not found: $REF"
    echo "FAILED missing_reference $(date)" >> "$STATUS_FILE"
    exit 1
fi

if [[ ! -f "${REF}.fai" ]]; then
    echo "ERROR: Reference FASTA index not found: ${REF}.fai"
    echo "FAILED missing_reference_index $(date)" >> "$STATUS_FILE"
    exit 1
fi


# Find source Octopus VCFs

mapfile -t PATIENT_VCFS < <(
    find "$VCF_BASE" \
        -type f \
        -name "${PATIENT}_*vsBLOOD_somatic.vcf.hg38_multianno.vcf" \
        | sort
)

if (( ${#PATIENT_VCFS[@]} == 0 )); then
    echo "ERROR: No source VCFs found for ${PATIENT}"
    echo "FAILED no_vcfs $(date)" >> "$STATUS_FILE"
    exit 1
fi

echo
echo "Found ${#PATIENT_VCFS[@]} source VCF(s):"
printf '  %s\n' "${PATIENT_VCFS[@]}"


# Find BAMs

mapfile -t PATIENT_BAMS < <(
    find "$BAM_BASE" \
        -maxdepth 2 \
        -type f \
        -name "XG-4050-${PATIENT}-*.md.bam" \
        | sort
)

if (( ${#PATIENT_BAMS[@]} == 0 )); then
    echo "ERROR: No BAMs found for ${PATIENT}"
    echo "FAILED no_bams $(date)" >> "$STATUS_FILE"
    exit 1
fi

echo
echo "Found ${#PATIENT_BAMS[@]} BAM(s):"
printf '  %s\n' "${PATIENT_BAMS[@]}"

# Check that BAM indices exist
for bam in "${PATIENT_BAMS[@]}"; do
    if [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" && ! -f "${bam}.csi" ]]; then
        echo "ERROR: BAM index not found for: $bam"
        echo "FAILED missing_bam_index $(date)" >> "$STATUS_FILE"
        exit 1
    fi
done

# Build union of Octopus SNVs

echo
echo "[$(date +%H:%M:%S)] Step 1: Building union of Octopus SNVs..."

TMP_DIR=$(mktemp -d "${OUT_DIR}/tmp.${PATIENT}.XXXXXX")
mkdir -p "${TMP_DIR}/sort"

NORMALIZED_VCFS=()

i=0
for v in "${PATIENT_VCFS[@]}"; do

    i=$((i + 1))

    out="${TMP_DIR}/source_${i}.vcf.gz"

    echo "  Processing: $(basename "$v")"

    bcftools view \
    -G \
    -v snps \
    -m2 \
    -M2 \
    -Ou \
    "$v" \
    | bcftools norm \
        -f "$REF" \
        --check-ref e \
        -Ou \
        - \
    | bcftools annotate \
        -x INFO \
        -Oz \
        -o "$out" \
        -


    bcftools index -t "$out"

    NORMALIZED_VCFS+=("$out")
done


# Combine all tissue-specific calls into one patient-level union.
# -d exact removes identical variants seen in more than one input VCF.

bcftools concat \
    -a \
    -d exact \
    "${NORMALIZED_VCFS[@]}" \
    -Ou \
    | bcftools sort \
        -T "${TMP_DIR}/sort" \
        -Ou \
    | bcftools norm \
        -d exact \
        -Oz \
        -o "$UNION_VCF"

bcftools index -t "$UNION_VCF"

n_union=$(bcftools view -H "$UNION_VCF" | wc -l)

if (( n_union == 0 )); then
    echo "ERROR: Union VCF contains no SNVs"
    echo "FAILED empty_union $(date)" >> "$STATUS_FILE"
    exit 1
fi

n_chr_union=$(
    bcftools query -f '%CHROM\n' "$UNION_VCF" \
        | sort -u \
        | wc -l
)

echo "  Union contains ${n_union} SNV loci across ${n_chr_union} chromosomes"


# Create allele target file for bcftools call -C alleles
# Format required by bcftools:
# CHROM  POS  REF,ALT

bcftools query \
    -f '%CHROM\t%POS\t%REF,%ALT\n' \
    "$UNION_VCF" \
    | bgzip -c \
    > "$ALLELES_TSV"

tabix -s1 -b2 -e2 "$ALLELES_TSV"


# Recover allele counts across all patient samples

echo
echo "[$(date +%H:%M:%S)] Step 2: Measuring union alleles in all BAMs..."

: > "$ERR_FILE"

if ! bcftools mpileup \
        -Ou \
        -f "$REF" \
        -R "$UNION_VCF" \
        -a FORMAT/AD,FORMAT/DP \
        -q 20 \
        -Q 20 \
        -d 5000 \
        -D \
        "${PATIENT_BAMS[@]}" \
        2> "$ERR_FILE" \
    | bcftools call \
        -m \
        -G - \
        -A \
        -C alleles \
        -T "$ALLELES_TSV" \
        -Oz \
        -o "$OUT_VCF" \
        2>> "$ERR_FILE"
then
    echo "ERROR: bcftools mpileup/call failed"
    tail -20 "$ERR_FILE"
    echo "FAILED forcecall $(date)" >> "$STATUS_FILE"
    exit 1
fi

bcftools index -t "$OUT_VCF"


# QC

n_out=$(bcftools view -H "$OUT_VCF" | wc -l)

n_chr_out=$(
    bcftools query -f '%CHROM\n' "$OUT_VCF" \
        | sort -u \
        | wc -l
)

echo "  Output contains ${n_out} loci across ${n_chr_out} chromosomes"

# Check that AD and DP exist
if ! bcftools query -l "$OUT_VCF" >/dev/null; then
    echo "ERROR: Output VCF cannot be read"
    echo "FAILED unreadable_output $(date)" >> "$STATUS_FILE"
    exit 1
fi

# The constrained output should contain the same target loci.
if (( n_out != n_union )); then
    echo "WARNING: Number of output loci (${n_out}) differs from target loci (${n_union})"
fi


# Final status

echo "DONE for ${PATIENT} at $(date)"
echo "Union SNV loci: ${n_union}"
echo "Output loci:    ${n_out}"
echo "Chromosomes:    ${n_chr_out}"

echo \
    "SUCCESS $(date) union_snvs=${n_union} output_records=${n_out} chrs=${n_chr_out}" \
    >> "$STATUS_FILE"
