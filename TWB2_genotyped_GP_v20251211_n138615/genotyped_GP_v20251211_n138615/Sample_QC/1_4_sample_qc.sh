# Combine Sex and HET QC results to remove samples failing either check
ARRAY=TWB2

awk "../Sex/"${ARRAY}"_S02_sex.problems.txt" "../HET/"${ARRAY}"_batch_remove_samples_het.txt" | unique | sort -u  > ${ARRAY}_QC_sex_het.txt

######################
# After GCR QC, do Sex check and het QC

plink2 \
--bfile GCR/"${ARRAY}"_S01_GCR \
--remove "${ARRAY}_QC_sex_het.txt" \
--make-bed \
--out "${ARRAY}_S03_QC
