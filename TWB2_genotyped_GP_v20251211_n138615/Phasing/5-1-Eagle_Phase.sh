
#!/bin/bash
set -Eeuo pipefail

# ==========================================
# TWB2 Eagle2 + Minimac4
# GRCh38 / 1000G EAS
# ==========================================

ARRAY="TWB2"

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

FAILED_FILE="${LOG_DIR}/failed_chr.txt"
OOM_FILE="${LOG_DIR}/oom_chr.txt"
SUCCESS_FILE="${LOG_DIR}/success_chr.txt"

: > "$FAILED_FILE"
: > "$OOM_FILE"
: > "$SUCCESS_FILE"

for chr in {1..22}; do

    echo "=================================="
    echo "Processing chr${chr}"
    echo "Start: $(date)"
    echo "=================================="

    INPUT="${INPUT_DIR}/${ARRAY}_chr${chr}.vcf.gz"

    REF_VCF="${REF_DIR}/chr${chr}.EAS.vcf.gz"
    REF_MSAV="${REF_DIR_MSAV}/chr${chr}.EAS.msav"

    # Eagle-compatible GRCh38 genetic map
    MAP="${MAP_DIR}/genetic_map_hg38_withX.txt.gz"

    PHASE_OUT="${OUT_DIR}/${ARRAY}_chr${chr}_phased"
    IMPUTE_OUT="${OUT_DIR}/${ARRAY}_chr${chr}_imputed"

    PHASE_VCF="${PHASE_OUT}.vcf.gz"

    PHASE_LOG="${LOG_DIR}/chr${chr}.phase.resource.log"
    IMPUTE_LOG="${LOG_DIR}/chr${chr}.impute.resource.log"

    # --------------------------------------
    # Input validation
    # --------------------------------------

    if [[ ! -s "$INPUT" ||
          ! -s "$REF_VCF" ||
          ! -s "$REF_MSAV" ||
          ! -s "$MAP" ]]; then

        echo "chr${chr}: Missing input file"
        echo "$chr" >> "$FAILED_FILE"
        continue
    fi

    # --------------------------------------
    # Step 1: Eagle Phasing
    # --------------------------------------

    if [[ -s "$PHASE_VCF" ]] && gzip -t "$PHASE_VCF" 2>/dev/null; then

        echo "chr${chr}: PHASED VCF EXISTS, SKIPPING EAGLE"

    else

        echo "chr${chr}: EAGLE PHASING"

        if ! /usr/bin/time -v \
            "$EAGLE" \
            --vcfRef "$REF_VCF" \
            --vcfTarget "$INPUT" \
            --geneticMapFile "$MAP" \
            --chrom "chr${chr}" \
            --allowRefAltSwap \
            --vcfOutFormat z \
            --numThreads "$THREADS" \
            --outPrefix "$PHASE_OUT" \
            > "$PHASE_LOG" 2>&1
        then

            echo "chr${chr}: PHASING FAILED"
            echo "$chr" >> "$FAILED_FILE"

            if grep -Eqi \
                "out of memory|bad_alloc|cannot allocate memory" \
                "$PHASE_LOG"; then
                echo "$chr:phasing" >> "$OOM_FILE"
            fi

            continue
        fi

    fi

    # --------------------------------------
    # Bcftools index the phased VCF
    # --------------------------------------

    bcftools index -f "$PHASE_VCF"
done

