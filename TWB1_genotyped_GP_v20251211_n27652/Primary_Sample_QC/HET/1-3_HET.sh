#!/bin/bash
# ==============================================================================
# Pipeline Step: TWB1 Heterozygosity Check
# Target Dataset: TWB1 cohort
# Inputs:
#   - TWB1_S01_GCR.{bed,bim,fam} (TWB1 dataset after GCR QC)
# Outputs:
#   - TWB1_S03_het_pending.{het,het.log} (Heterozygosity check results)
# Do:
#   1. Perform heterozygosity check using PLINK2
# ==============================================================================

WORKING_PREFIX=../GCR/TWB1_S01_GCR
ARRAY=TWB1

plink2 \
--bfile "${WORKING_PREFIX}" \
--autosome \
--het \
--out "${ARRAY}"_S03_het_pending