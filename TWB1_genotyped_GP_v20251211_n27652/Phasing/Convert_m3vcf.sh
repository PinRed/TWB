#!/bin/bash
set -Eeuo pipefail
INPUT="1kgp_v5/EAS_filtered"
OUTPUT="1kgp_v5/Minimac_ReHap"
mkdir -p ${OUTPUT}

for chr in {1..22}; do
    echo "=================================="
    echo "Processing chr${chr}"
    echo "Start: $(date)"
    echo "=================================="

    # Convert VCF to m3vcf for Minimac4
    echo "chr${chr}: Converting VCF to m3vcf for Minimac4"

    minimac4 \
        --compress-reference \
        ${INPUT}/chr${chr}.EAS.vcf.gz \
        -o ${OUTPUT}/chr${chr}.EAS.msav
done

    