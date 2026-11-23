#!/bin/bash
set -euo pipefail

Input_dir="1kgp_v5"
Output_dir="1kgp_v5/EAS_filtered"

mkdir -p "${Output_dir}"

INFO="/SATA2TB/Taiwan_Biobank/1kGP_SGDP_prep/metadata/igsr_samples.tsv"

# Step 1: Extract EAS sample IDs
awk 'NR>1 && $6=="EAS" {print $1}' "$INFO" \
    > "${Output_dir}/EAS.id"

awk 'NR==FNR{a[$0];next} $0 in a' \
    "${Output_dir}/EAS.id" \
    "${Output_dir}/keep_unrelated_samples.txt" \
    > "${Output_dir}/EAS_unrelated.id"

echo "EAS samples:"
wc -l "${Output_dir}/EAS_unrelated.id"

# Step 2: Process chr1–22
for chr in {1..22}; do

    echo "Processing chromosome ${chr}..."

    input="${Input_dir}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
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
        -e 'INFO/SVTYPE!="."' \
        -Ou |

    bcftools view \
        -S "${Output_dir}/EAS_unrelated.id" \
        -Ou |

    bcftools +fill-tags \
        -Ou -- -t AC,AN,AF |

    bcftools view \
        -i 'INFO/AN>=4 && (SUM(INFO/AC)>=2 && INFO/AN-SUM(INFO/AC)>=2)' \
        -Oz -o "$output"

    rm "${Output_dir}/contig.${chr}.txt"

    else

        bcftools view \
            -e 'INFO/SVTYPE!="."' \
            -Ou "$input" |

        bcftools view \
            -S "${Output_dir}/EAS_unrelated.id" \
            -Ou |

        bcftools +fill-tags \
            -Ou -- -t AC,AN,AF |

        bcftools view \
            -i 'INFO/AN>=4 && (SUM(INFO/AC)>=2 && INFO/AN-SUM(INFO/AC)>=2)' \
            -Oz -o "$output"

    fi

    # Remove non-unique IDs (set ID to ".")
    bcftools query -f '%ID\n' "$output" |
        tr ';' '\n' |
        awk '$0!="." && $0!=""' |
        sort |
        uniq -d > "${Output_dir}/duplicate.${chr}.txt"

    bcftools view -h "$output" > "${output}.tmp.vcf"

    bcftools view -H "$output" |
        awk -F '\t' '
            BEGIN {OFS="\t"}

            FILENAME==ARGV[1] {
                dup[$1]=1
                next
            }

            {
                n=split($3, ids, ";")

                for (i=1; i<=n; i++) {
                    if (ids[i] in dup) {
                        $3="."
                        break
                    }
                }

                print
            }
        ' "${Output_dir}/duplicate.${chr}.txt" - \
        >> "${output}.tmp.vcf"

    bcftools view \
        -Oz -o "${output}.filtered.vcf.gz" \
        "${output}.tmp.vcf"

    mv "${output}.filtered.vcf.gz" "$output"

    rm "${output}.tmp.vcf"
    #rm "${Output_dir}/duplicate.${chr}.txt"

    # Index
    bcftools index -f -t "$output"

    echo "Chromosome ${chr} completed."

done

echo "All chromosomes completed."