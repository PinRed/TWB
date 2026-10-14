## Linux
PREFIX="../../TWB2.hg38.GP.20251211"
# TWB 2 Batch sample missing rate 
for BATCH in {0001..0154}; do
    KEEP="../../Sample_Batch/TWB2_B${BATCH}_samples.keep"

    if [[ ! -s "$KEEP" ]]; then
        echo "Skipping batch ${BATCH}: ${KEEP} missing or empty"
        continue
    fi

    plink2 \
        --bfile "$PREFIX" \
        --keep "$KEEP" \
        --missing \
        --out "Missing/TWB2_Batch${BATCH}_missing"
done

# TWB_Sample_QC.R
# output TWB2_S01_GCR_exclude.txt

plink2 --bfile "$PREFIX" --remove TWB2_S01_GCR_exclude.txt --make-bed --out TWB2_S01_GCR