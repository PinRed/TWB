#!/bin/bash
# ==============================================================================
# Pipeline Step: TWB2 Sex & HET QC
# Target Dataset: TWB2 cohort
# Inputs:
#   - TWB2_S01_GCR.{bed,bim,fam} (TWB2 dataset after GCR QC)
#   - TWB2_S02_sex.{sexcheck,summary.tsv,problems.txt} (Sex check results and summary)
#   - TWB2_S03_het_pending.{het,het.log} (Heterozygosity check results)
# Outputs:
#   - TWB2_S03_QC.{bed,bim,fam} (TWB2 dataset after removing samples failing Sex or HET QC)
# Do:
#   1. Combine Sex and HET QC results to remove samples failing either check
#   2. Create a new dataset after removing the failed samples
# ==============================================================================
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
