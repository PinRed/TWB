## Linux
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