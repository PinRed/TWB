WORKING_PREFIX=../GCR/TWB1_S01_GCR
ARRAY=TWB1

plink2 \
--bfile "${WORKING_PREFIX}" \
--autosome \
--het \
--out "${ARRAY}"_S03_het_pending