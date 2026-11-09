#!/bin/bash
set -euo pipefail

Input_dir="1kgp_v5"
Output_dir="1kgp_v5/EAS_filtered"

mkdir -p "${Output_dir}"

INFO="${Input_dir}/integrated_call_samples_v3.20130502.ALL.panel"

# Step 1: Extract EAS sample IDs
awk 'NR>1 && $3=="EAS" {print $1}' "$INFO" \
    > "${Output_dir}/EAS.id"

echo "EAS samples:"
wc -l "${Output_dir}/EAS.id"

# Step 2: Process chr1–22
for chr in {1..22}; do

    echo "Processing chromosome ${chr}..."

    input="${Input_dir}/chr${chr}.1kg.phase3.v5a.vcf.gz"
    output="${Output_dir}/chr${chr}.EAS.vcf.gz"

    # Check chromosome naming and contig header
    if ! bcftools view -h "$input" |
        grep -q "^##contig=<ID=${chr}[,>]"; then

        printf '##contig=<ID=%s>\n' "$chr" \
            > "${Output_dir}/contig.${chr}.txt"

        bcftools annotate \
            -h "${Output_dir}/contig.${chr}.txt" \
            -Ou "$input" |

        bcftools view \
            -S "${Output_dir}/EAS.id" \
            -Ou |

        bcftools +fill-tags \
            -Ou -- -t AC,AN,AF |

        bcftools view \
            -e 'INFO/AN=0 || INFO/AC=0 || INFO/AC=INFO/AN' \
            -Oz -o "$output"

        rm "${Output_dir}/contig.${chr}.txt"

    else

        bcftools view \
            -S "${Output_dir}/EAS.id" \
            -Ou "$input" |

        bcftools +fill-tags \
            -Ou -- -t AC,AN,AF |

        bcftools view \
            -e 'INFO/AN=0 || INFO/AC=0 || INFO/AC=INFO/AN' \
            -Oz -o "$output"

    fi

    # Index
    bcftools index -t "$output"

    echo "Chromosome ${chr} completed."

done

echo "All chromosomes completed."