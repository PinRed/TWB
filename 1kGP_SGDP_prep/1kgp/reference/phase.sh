wget -c \
'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz'

wget -c \
'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz.tbi'

bcftools query \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
  ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38.27022019.sites.vcf.gz \
  > 1000G_Phase3_GRCh38.tsv


awk -F '\t' '
BEGIN {OFS="\t"}
{
    if (length($4)==1 && length($5)==1)
        snp++
    else
        indel++
    if ($3==".")
        no_id++
}
END {
    print "SNP:", snp
    print "INDEL/MNV:", indel
    print "No ID:", no_id
}' 1000G_Phase3_GRCh38.tsv

awk -F '\t' '
length($4)==1 && length($5)==1 &&
$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {
    print
}' 1000G_Phase3_GRCh38.tsv \
> 1000G_Phase3_GRCh38_SNP.tsv