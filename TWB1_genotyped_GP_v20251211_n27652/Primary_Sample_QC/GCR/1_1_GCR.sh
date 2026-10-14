#!/bin/bash
# ==============================================================================
# Pipeline Step: TWB2 Batch Sample Missing Rate Calculation
# Target Dataset: TWB1 cohort 
# Inputs:
#   - TWB1.hg19.GP.20251211.{bed,bim,fam} (TWB1 genotyped dataset)
#   - TWB1_B{1..14}_samples.keep (List of samples for each batch)
# Outputs:
#   - Missing/TWB1_Batch{1..14}_missing.{imiss,imiss.log} (Missing rate for each batch)
#   - TWB1_S01_GCR.{bed,bim,fam} (TWB1 dataset after removing samples failing GCR QC)
# Do:
#   1. Calculate missing rate for each batch of TWB1 samples
#   2. Remove samples failing GCR QC and create a new dataset
# ==============================================================================

PREFIX="../../TWB1.hg19.GP.20251211"
# TWB 1 Batch sample missing rate 
for BATCH in {1..14}; do
    KEEP="../../Sample_Batch/TWB1_B${BATCH}_samples.keep"

    if [[ ! -s "$KEEP" ]]; then
        echo "Skipping batch ${BATCH}: ${KEEP} missing or empty"
        continue
    fi

    plink2 \
        --bfile "$PREFIX" \
        --keep "$KEEP" \
        --missing \
        --out "Missing/TWB1_Batch${BATCH}_missing"
done

# TWB_Sample_QC.R
# output TWB1_S01_GCR_exclude.txt

plink2 --bfile "$PREFIX" --remove TWB1_S01_GCR_exclude.txt --make-bed --out TWB1_S01_GCR