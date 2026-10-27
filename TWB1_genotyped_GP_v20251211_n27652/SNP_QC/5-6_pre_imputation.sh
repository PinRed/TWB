#!/bin/bash
set -euo pipefail

# =================== Paths ========================
ARRAY="TWB1"
INPUT_DIR="SNP_QC_FILE"
HARM_DIR="Harmonization"
BFILE="${INPUT_DIR}/${ARRAY}_F_SNPQC05_clean"

# hg19 reference FASTA（請修改為實際路徑）
REF_FA="GRCh37_genome.fa"

# SNP QC lists
EXCLUDE="${HARM_DIR}/approved_exclude.ids"
FLIP="${HARM_DIR}/approved_flip.ids"

# Output directory
OUT="${HARM_DIR}"
OUT_FILE_NAME="${ARRAY}_F_SNPQC06"
mkdir -p "$OUT"

BASE="${OUT}/BASE"
ELIGIBLE="${OUT}/${OUT_FILE_NAME}_eligible"
FLIPPED="${OUT}/${OUT_FILE_NAME}_flipped"
REFASSIGNED="${OUT}/${OUT_FILE_NAME}_refassigned"

# =================== Check files ========================

for ext in bed bim fam; do
    [[ -f "${BFILE}.${ext}" ]] || {
        echo "ERROR: Missing ${BFILE}.${ext}"
        exit 1
    }
done

[[ -f "$REF_FA" ]] || {
    echo "ERROR: Missing FASTA: $REF_FA"
    exit 1
}

[[ -f "$EXCLUDE" ]] || {
    echo "ERROR: Missing exclude list"
    exit 1
}

[[ -f "$FLIP" ]] || {
    echo "ERROR: Missing flip list"
    exit 1
}

# =================== Step 1: BED -> PGEN ========================

echo "Step 1: Convert BED to PGEN"

plink2 \
    --bfile "$BFILE" \
    --make-pgen \
    --out "$BASE"

# =================== Step 2: Exclude SNPs ========================

echo "Step 2: Exclude SNPs"

plink2 \
    --pfile "$BASE" \
    --exclude "$EXCLUDE" \
    --make-pgen \
    --out "$ELIGIBLE"

# =================== Step 3: Strand flip ========================

echo "Step 3: Strand flip"

if [[ -s "$FLIP" ]]; then

    plink2 \
        --pfile "$ELIGIBLE" \
        --flip "$FLIP" \
        --make-pgen \
        --out "$FLIPPED"

else

    echo "No approved flip SNPs. Skipping flip."

    cp "${ELIGIBLE}.pgen" "${FLIPPED}.pgen"
    cp "${ELIGIBLE}.pvar" "${FLIPPED}.pvar"
    cp "${ELIGIBLE}.psam" "${FLIPPED}.psam"

fi

# =================== Step 4: Reference assignment ========================

echo "Step 4: Assign REF allele from hg19 FASTA"

plink2 \
    --pfile "$FLIPPED" \
    --fa "$REF_FA" \
    --ref-from-fa \
    --make-pgen \
    --out "$REFASSIGNED"

# =================== Step 5: QC Summary ========================

echo ""
echo "========== Harmonization Summary =========="

for prefix in "$BASE" "$ELIGIBLE" "$FLIPPED" "$REFASSIGNED"; do

    N=$(awk '!/^#/ {n++} END {print n+0}' "${prefix}.pvar")

    echo "$(basename "$prefix"): $N variants"

done

echo ""
echo "Excluded SNP IDs: $(wc -l < "$EXCLUDE")"
echo "Flip SNP IDs:     $(wc -l < "$FLIP")"

echo ""
echo "Final output:"
echo "${REFASSIGNED}.pgen"
echo "${REFASSIGNED}.pvar"
echo "${REFASSIGNED}.psam"

echo "========== Completed =========="

rm -rf "$BASE".* "$ELIGIBLE".* "$FLIPPED".*

# ============================================================
# Step 6 — Build validated chromosome VCF
# ============================================================

echo ""
echo "========== Step 6: Build validated VCF =========="

VCF_DIR="${HARM_DIR}/Validated_VCF"
mkdir -p "$VCF_DIR"

# TWB1: GRCh37, chromosome 1-22
# Pilot: CHROMS=(20)
# All chromosomes:
CHROMS=($(seq 1 22))

# Required tools
for tool in plink2 bcftools sha256sum awk sort uniq comm; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: Missing tool: $tool"
        exit 1
    }
done

# ============================================================
# Step 7 — Sample ID QC
# ============================================================

echo "========== Step 7: Sample ID QC =========="

SAMPLE_DIR="${VCF_DIR}/sample_qc"
mkdir -p "$SAMPLE_DIR"

awk '
BEGIN {FS="[ \t]+"}
NR==1 {
    for (i=1; i<=NF; i++) {
        col=$i
        sub(/^#/, "", col)
        if (col=="IID") iid_col=i
    }
    if (!iid_col) exit 2
    next
}
!/^#/ && NF>0 {
    print $iid_col
}
' "${REFASSIGNED}.psam" \
    > "${SAMPLE_DIR}/IID.txt" || {
        echo "ERROR: Cannot identify IID in PSAM"
        exit 1
    }

SAMPLE_N=$(wc -l < "${SAMPLE_DIR}/IID.txt")

sort "${SAMPLE_DIR}/IID.txt" |
    uniq -d > "${SAMPLE_DIR}/duplicate_IID.txt"

if [[ -s "${SAMPLE_DIR}/duplicate_IID.txt" ]]; then
    echo "ERROR: Duplicate IID detected"
    exit 1
fi

sort -u "${SAMPLE_DIR}/IID.txt" \
    > "${SAMPLE_DIR}/IID.sorted.txt"

IID_SHA256=$(sha256sum "${SAMPLE_DIR}/IID.sorted.txt" | awk '{print $1}')

echo "Sample N: ${SAMPLE_N}"
echo "IID SHA256: ${IID_SHA256}"

# ============================================================
# Step 8 — Manifest
# ============================================================

MANIFEST="${VCF_DIR}/${ARRAY}_VCF_manifest.tsv"

printf "ARRAY\tBUILD\tCHR\tSAMPLE_N\tIID_SHA256\tVARIANT_N\tDUPLICATE_POS_N\tREF_STATUS\tVCF_SHA256\tVCF\n" \
    > "$MANIFEST"

# ============================================================
# Step 9 — Export, sort, index, validate each chromosome
# ============================================================

for CHR in "${CHROMS[@]}"; do

    echo ""
    echo "========================================"
    echo "Processing ${ARRAY} chromosome ${CHR}"
    echo "========================================"

    CHR_DIR="${VCF_DIR}/chr${CHR}"
    mkdir -p "$CHR_DIR"

    PREFIX="${CHR_DIR}/${ARRAY}_chr${CHR}"
    RAW="${PREFIX}_raw.vcf.gz"
    FINAL="${PREFIX}.vcf.gz"

    # --------------------------------------------------------
    # 9.1 Check chromosome exists in GRCh37 FASTA
    # --------------------------------------------------------

    awk -v c="$CHR" '$1==c {found=1}
        END {exit !found}' "${REF_FA}.fai" || {
        echo "ERROR: Chromosome ${CHR} absent from FASTA"
        exit 1
    }

    # --------------------------------------------------------
    # 9.2 Export chromosome VCF
    # --------------------------------------------------------

    echo "[1/7] Export chromosome ${CHR}"

    plink2 \
        --pfile "$REFASSIGNED" \
        --chr "$CHR" \
        --output-chr 26 \
        --export vcf bgz id-paste=iid\
        --out "${PREFIX}_raw"

    [[ -s "$RAW" ]] || {
        echo "ERROR: Raw VCF missing"
        exit 1
    }

    # --------------------------------------------------------
    # 9.3 Sort VCF
    # --------------------------------------------------------

    echo "[2/7] Sort VCF"

    bcftools sort "$RAW" \
        -Oz \
        -o "$FINAL"

    # --------------------------------------------------------
    # 9.4 Index VCF
    # --------------------------------------------------------

    echo "[3/7] Index VCF"

    bcftools index \
        -f -t "$FINAL"

    [[ -s "${FINAL}.tbi" ]] || {
        echo "ERROR: TBI index missing"
        exit 1
    }

    # --------------------------------------------------------
    # 9.5 Sample QC
    # --------------------------------------------------------

    echo "[4/7] Sample ID validation"

    bcftools query -l "$FINAL" \
        > "${PREFIX}.samples.txt"

    VCF_SAMPLE_N=$(wc -l < "${PREFIX}.samples.txt")

    [[ "$VCF_SAMPLE_N" -eq "$SAMPLE_N" ]] || {
        echo "ERROR: Sample count mismatch"
        exit 1
    }

    sort -u "${PREFIX}.samples.txt" \
        > "${PREFIX}.samples.sorted.txt"

    comm -3 \
        "${SAMPLE_DIR}/IID.sorted.txt" \
        "${PREFIX}.samples.sorted.txt" \
        > "${PREFIX}.sample_diff.txt"

    [[ ! -s "${PREFIX}.sample_diff.txt" ]] || {
        echo "ERROR: Sample ID mismatch"
        exit 1
    }

    # --------------------------------------------------------
    # 9.6 Variant count and duplicate POS
    # --------------------------------------------------------

    echo "[5/7] Variant QC"

    VARIANT_N=$(bcftools view -H "$FINAL" | wc -l)

    [[ "$VARIANT_N" -gt 0 ]] || {
        echo "ERROR: No variants on chromosome ${CHR}"
        exit 1
    }

    bcftools query \
        -f '%CHROM\t%POS\n' \
        "$FINAL" |
        sort |
        uniq -d \
        > "${PREFIX}.duplicate_pos.txt"

    DUP_POS_N=$(wc -l < "${PREFIX}.duplicate_pos.txt")

    if [[ "$DUP_POS_N" -gt 0 ]]; then
        echo "ERROR: Duplicate POS on chr${CHR}: ${DUP_POS_N}"
        exit 1
    fi

    # --------------------------------------------------------
    # 9.7 REF validation against GRCh37 FASTA
    # --------------------------------------------------------

    echo "[6/7] REF validation"

    bcftools norm \
        --check-ref e \
        -f "$REF_FA" \
        "$FINAL" \
        -Ou \
        -o /dev/null \
        2> "${PREFIX}.ref_check.log" || {
        echo "ERROR: REF validation failed on chr${CHR}"
        exit 1
    }

    REF_STATUS="PASS"

    # --------------------------------------------------------
    # 9.8 SHA256 and manifest
    # --------------------------------------------------------

    echo "[7/7] SHA256 and manifest"

    VCF_SHA256=$(sha256sum "$FINAL" | awk '{print $1}')

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$ARRAY" \
        "GRCh37" \
        "$CHR" \
        "$SAMPLE_N" \
        "$IID_SHA256" \
        "$VARIANT_N" \
        "$DUP_POS_N" \
        "$REF_STATUS" \
        "$VCF_SHA256" \
        "$FINAL" \
        >> "$MANIFEST"

    # Remove raw VCF only after all checks pass
    rm -f "$RAW"

    echo "PASS: chr${CHR}"
    echo "Variants: ${VARIANT_N}"
    echo "Samples: ${SAMPLE_N}"
    echo "REF: ${REF_STATUS}"

done

# ============================================================
# Step 10 — Final summary
# ============================================================

echo ""
echo "========================================"
echo "VCF VALIDATION COMPLETED"
echo "========================================"

echo "Array: ${ARRAY}"
echo "Build: GRCh37"
echo "Chromosomes: ${CHROMS[*]}"
echo "Sample N: ${SAMPLE_N}"
echo "Manifest: ${MANIFEST}"
echo "Output directory: ${VCF_DIR}"

echo ""
echo "All requested chromosomes PASSED."

echo "Move the final VCFs "

mkdir -p "${HARM_DIR}/VCF_Files"
OUTPUT_VCF_DIR="${HARM_DIR}/VCF_Files"

for CHR in "${CHROMS[@]}"; do
    echo "Move  ${VCF_DIR}/chr${CHR}/${ARRAY}_chr${CHR}.vcf.gz"
    mv "${VCF_DIR}/chr${CHR}/${ARRAY}_chr${CHR}.vcf.gz" "${OUTPUT_VCF_DIR}/"
    mv "${VCF_DIR}/chr${CHR}/${ARRAY}_chr${CHR}.vcf.gz.tbi" "${OUTPUT_VCF_DIR}/"
done

echo "Removing temporary directories"
rm -rf "${VCF_DIR}"