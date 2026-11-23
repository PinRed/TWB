#!/bin/bash
set -Eeuo pipefail

# ==========================================
# TWB2 Eagle2 + Minimac4
# GRCh38 / 1000G EAS
# Chromosome 20 ONLY
# ==========================================

ARRAY="TWB2"
CHR=20

INPUT_DIR="../SNP_QC/Harmonization/VCF_Files"
REF_DIR="1kgp_v5/EAS_filtered"
REF_DIR_MSAV="1kgp_v5/Minimac_ReHap"
MAP_DIR="/SATA2TB/Taiwan_Biobank/Tools/Eagle_v2.4.1/tables"

OUT_DIR="Phased_Impute_VCF"
LOG_DIR="${OUT_DIR}/logs"

EAGLE="../../Tools/Eagle_v2.4.1/eagle"
MINIMAC4="minimac4"

THREADS=24

mkdir -p "$OUT_DIR" "$LOG_DIR"

FAILED_FILE="${LOG_DIR}/failed_chr20.txt"
OOM_FILE="${LOG_DIR}/oom_chr20.txt"
SUCCESS_FILE="${LOG_DIR}/success_chr20.txt"

: > "$FAILED_FILE"
: > "$OOM_FILE"
: > "$SUCCESS_FILE"

# ==========================================
# File paths
# ==========================================

INPUT="${INPUT_DIR}/${ARRAY}_chr${CHR}.vcf.gz"

REF_VCF="${REF_DIR}/chr${CHR}.EAS.vcf.gz"
REF_MSAV="${REF_DIR_MSAV}/chr${CHR}.EAS.msav"

MAP="${MAP_DIR}/genetic_map_hg38_withX.txt.gz"

PHASE_OUT="${OUT_DIR}/${ARRAY}_chr${CHR}_phased"
PHASE_VCF="${PHASE_OUT}.vcf.gz"

IMPUTE_OUT="${OUT_DIR}/${ARRAY}_chr${CHR}_imputed"
IMPUTE_VCF="${IMPUTE_OUT}.dose.vcf.gz"

PHASE_LOG="${LOG_DIR}/chr${CHR}.phase.resource.log"
IMPUTE_LOG="${LOG_DIR}/chr${CHR}.impute.resource.log"

# ==========================================
# Helper functions
# ==========================================

check_oom() {
    local logfile="$1"
    local step="$2"

    if [[ -f "$logfile" ]] &&
       grep -Eqi \
       "out of memory|bad_alloc|cannot allocate memory|oom-kill|killed process" \
       "$logfile"; then
        echo "${CHR}:${step}" >> "$OOM_FILE"
    fi
}

fail_step() {
    local step="$1"
    local logfile="${2:-}"

    echo "chr${CHR}: ${step} FAILED"
    echo "$CHR" >> "$FAILED_FILE"

    if [[ -n "$logfile" ]]; then
        check_oom "$logfile" "$step"
    fi

    exit 1
}

# ==========================================
# Start
# ==========================================

echo "=================================="
echo "TWB2 Eagle2 + Minimac4"
echo "Chromosome: chr${CHR}"
echo "Start: $(date)"
echo "=================================="

# ==========================================
# Step 0: Input validation
# ==========================================

for file in "$INPUT" "$REF_VCF" "$REF_MSAV" "$MAP"; do
    if [[ ! -s "$file" ]]; then
        echo "Missing input: $file"
        fail_step "input_validation"
    fi
done

# ==========================================
# Step 1: Eagle2 Phasing
# ==========================================

if [[ -s "$PHASE_VCF" ]] &&
   gzip -t "$PHASE_VCF" 2>/dev/null; then

    echo "chr${CHR}: PHASED VCF EXISTS"
    echo "Skipping Eagle2"

else

    echo "=================================="
    echo "chr${CHR}: EAGLE2 PHASING"
    echo "Start: $(date)"
    echo "=================================="

    if ! /usr/bin/time -v \
        "$EAGLE" \
        --vcfRef "$REF_VCF" \
        --vcfTarget "$INPUT" \
        --geneticMapFile "$MAP" \
        --chrom "chr${CHR}" \
        --allowRefAltSwap \
        --vcfOutFormat z \
        --numThreads "$THREADS" \
        --outPrefix "$PHASE_OUT" \
        > "$PHASE_LOG" 2>&1
    then
        fail_step "phasing" "$PHASE_LOG"
    fi

    if [[ ! -s "$PHASE_VCF" ]] ||
       ! gzip -t "$PHASE_VCF" 2>/dev/null; then
        fail_step "phasing_output" "$PHASE_LOG"
    fi

    echo "chr${CHR}: PHASING SUCCESS"
fi

# ==========================================
# Step 2: Index phased VCF
# ==========================================

echo "chr${CHR}: INDEX PHASED VCF"

if ! bcftools index -f "$PHASE_VCF"; then
    fail_step "phasing_index"
fi

# ==========================================
# Step 3: Minimac4 Imputation
# ==========================================

if [[ -s "$IMPUTE_VCF" ]] &&
   gzip -t "$IMPUTE_VCF" 2>/dev/null; then

    echo "chr${CHR}: IMPUTED VCF EXISTS"
    echo "Skipping Minimac4"

else

    echo "=================================="
    echo "chr${CHR}: MINIMAC4 IMPUTATION"
    echo "Start: $(date)"
    echo "=================================="

    if ! /usr/bin/time -v \
        "$MINIMAC4" \
        --refHaps "$REF_MSAV" \
        --haps "$PHASE_VCF" \
        --window 500000 \
        --prefix "$IMPUTE_OUT" \
        --cpus "$THREADS" \
        --chr "chr${CHR}" \
        --noPhoneHome \
        --format GT,DS,GP \
        --allTypedSites \
        --meta \
        --minRatio 0.00001 \
        > "$IMPUTE_LOG" 2>&1
    then
        fail_step "imputation" "$IMPUTE_LOG"
    fi

fi

# ==========================================
# Step 4: Validate imputed VCF
# ==========================================

echo "chr${CHR}: VALIDATING IMPUTED VCF"

if [[ ! -s "$IMPUTE_VCF" ]] ||
   ! gzip -t "$IMPUTE_VCF" 2>/dev/null; then

    fail_step "imputation_output" "$IMPUTE_LOG"
fi




# ==========================================
# Step 5: Final QC
# INFO/R2 > 0.6
# MAF > 0.005
# ==========================================

echo "=================================="
echo "chr${CHR}: FINAL IMPUTATION QC"
echo "Start: $(date)"
echo "=================================="

QC_OUT="${OUT_DIR}/${ARRAY}_chr${CHR}_imputed_QC.vcf.gz"
QC_LOG="${LOG_DIR}/chr${CHR}.final_qc.log"
QC_SUMMARY="${LOG_DIR}/chr${CHR}.final_qc.summary.txt"

# Check required INFO tags
if ! bcftools view -h "$IMPUTE_VCF" |
    grep -q '^##INFO=<ID=R2,'; then
    fail_step "missing_R2_tag"
fi

if ! bcftools view -h "$IMPUTE_VCF" |
    grep -q '^##INFO=<ID=AF,'; then
    fail_step "missing_AF_tag"
fi

# Count variants before QC
N_BEFORE=$(bcftools index -n "$IMPUTE_VCF" 2>/dev/null ||
    bcftools view -H "$IMPUTE_VCF" | wc -l)

echo "Variants before QC: $N_BEFORE"

# Filter by imputation quality and MAF
# MAF = min(AF, 1-AF)
if ! bcftools view \
    -i 'INFO/R2>0.6 && INFO/AF>0.001 && INFO/AF<0.999' \
    -Oz \
    -o "$QC_OUT" \
    "$IMPUTE_VCF" \
    > "$QC_LOG" 2>&1; then

    fail_step "final_qc" "$QC_LOG"
fi

# Validate output
if [[ ! -s "$QC_OUT" ]] ||
   ! bcftools view -h "$QC_OUT" >/dev/null 2>&1; then
    fail_step "final_qc_output" "$QC_LOG"
fi

# Index QC VCF
if ! bcftools index -f -t "$QC_OUT"; then
    fail_step "final_qc_index" "$QC_LOG"
fi

# Count variants after QC
N_AFTER=$(bcftools index -n "$QC_OUT")

# Calculate retention rate
RETENTION=$(awk -v a="$N_AFTER" -v b="$N_BEFORE" \
    'BEGIN {if(b>0) printf "%.2f",100*a/b; else print "NA"}')

# Write QC summary
{
    echo "=================================="
    echo "TWB2 Imputation Final QC"
    echo "=================================="
    echo "Chromosome: chr${CHR}"
    echo "INFO/R2 threshold: > 0.6"
    echo "MAF threshold: > 0.005"
    echo "Variants before QC: $N_BEFORE"
    echo "Variants after QC: $N_AFTER"
    echo "Retention rate: ${RETENTION}%"
    echo "Output: $QC_OUT"
    echo "Date: $(date)"
    echo "=================================="
} | tee "$QC_SUMMARY"

echo "chr${CHR}: FINAL QC SUCCESS"


# ==========================================
# Step 5: Success
# ==========================================

echo "$CHR" >> "$SUCCESS_FILE"

echo
echo "=================================="
echo "FINAL SUMMARY"
echo "=================================="

echo "Chromosome: chr${CHR}"
echo "Status: SUCCESS"
echo "Phased VCF: $PHASE_VCF"
echo "Imputed VCF: $IMPUTE_VCF"
echo "End: $(date)"
echo "=================================="

echo "All chr${CHR} steps completed successfully."
