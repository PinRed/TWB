# Genomic Data Sources

This repository uses genomic data from the following datasets:

- Taiwan Biobank (TWB)
- 1000 Genomes Project (1KGP)
- Simons Genome Diversity Project (SGDP)

## 1. Taiwan Biobank (TWB)

- **Data source:** [Taiwan Biobank](https://www.biobank.org.tw/)
- **TWB1**
  - Reference genome: `GRCh37`
  - Reference publication: [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2666979X2200146X)
- **TWB2**
  - Reference genome: `GRCh38`
  - Reference publication: [Nature](https://www.nature.com/articles/s41586-025-09680-x)

## 2. 1000 Genomes Project (1KGP)

- **GRCh37**
  - Data release: `20130502`
  - Download: [1000 Genomes GRCh37](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)
- **GRCh38**
  - Data collection: `1000G_2504_high_coverage`
  - Phased dataset: `20201028_3202_phased`
  - Download: [1000 Genomes GRCh38](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)

## 3. Simons Genome Diversity Project (SGDP)

- **GRCh37**
  - Data format: `VCF`
  - Download: [SGDP VCF Variants](https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/vcf_variants/)
 
# TWB Genotype Quality Control Pipeline

## Quality Control Summary

Summary of genotype quality control (QC) procedures for Taiwan Biobank Phase 1 (TWB1) and Phase 2 (TWB2).

| Step | TWB1 | TWB2 | QC Criteria / Notes | Output File |
|---|---|---|---|---|
| **Original** | Samples: 27,652<br>SNPs: 646,973 | Samples: 138,615<br>SNPs: 686,370 | Raw genotype data | - |
| **Primary Sample QC** | Samples: 27,323<br>SNPs: 646,973 | Samples: 136,984<br>SNPs: 646,973 | GCR, SEX, HET | - |
| **Primary SNP QC** | Samples: 27,323<br>SNPs: 611,991 | Samples: 138,615<br>SNPs: 571,914 | Structural SNP filtering:<br>MAF = 0<br>Missingness < 98% | `TWBX_S03_all_QC_clean` |
| **EAS Kinship Filtering** | Samples: 25,206<br>SNPs: 611,898 | Samples: 136,984<br>SNPs: 539,122 | Remove batch-associated SNPs<br>Remove related samples (Kinship > 0.0884) | `TWBX_S05` |
| **Meta PCA** | Samples: 25,206<br>SNPs: 62,944 | Samples: 136,984<br>SNPs: 47,514 | Combined datasets: 1KGP + SGDP + TWB | `TWBX_S05_pca` |
| **EAS PCAiR** | Samples: 25,205<br>SNPs: 147,313 | Samples: 136,979<br>SNPs: 142,833 | Select EAS samples<br>Remove batch-associated SNPs<br>MAF < 0.05<br>Missingness < 99%<br>LD pruning: 5000 kb, r² = 0.3 | `TWBX__EAS_pcair_kin` |
| **Final SNP QC** | Samples: 25,205<br>SNPs: 586,508 | Samples: 136,979<br>SNPs: 479,138 | Start from Primary SNP QC<br>Select EAS samples<br>SNP missingness < 98%<br>MAF = 0<br>HWE < 1e-6<br>Batch-effect SNP filtering | `TWBX_F_SNPQC05_clean.bim` |
