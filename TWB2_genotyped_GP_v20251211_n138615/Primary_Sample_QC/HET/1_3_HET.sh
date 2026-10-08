WORKING_PREFIX=../GCR/TWB2_S01_GCR
ARRAY=TWB2

plink2 \
--bfile "${WORKING_PREFIX}" \
--autosome \
--het \
--out "${ARRAY}"_S03_het_pending