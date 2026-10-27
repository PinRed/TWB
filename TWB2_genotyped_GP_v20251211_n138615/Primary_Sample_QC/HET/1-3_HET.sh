#!/bin/bash
# ==============================================================================
# Pipeline Step: TWB2 Heterozygosity Check
# Target Dataset: TWB2 cohort
# Inputs:
#   - TWB2_S01_GCR.{bed,bim,fam} (TWB2 dataset after GCR QC)
# Outputs:
#   - TWB2_S03_het_pending.{het,het.log} (Heterozygosity check results)
# Do:
#   1. Perform heterozygosity check using PLINK2
# ==============================================================================
WORKING_PREFIX=../GCR/TWB2_S01_GCR
ARRAY=TWB2

plink2 \
--bfile "${WORKING_PREFIX}" \
--autosome \
--het \
--out "${ARRAY}"_S03_het_pending