wget -c \
'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz'

wget -c \
'https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz.tbi'

bcftools query \
  -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' \
  ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz \
  > 1000G_Phase3_GRCh37.tsv


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
}' 1000G_Phase3_GRCh37.tsv

awk -F '\t' '
length($4)==1 && length($5)==1 &&
$4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/ {
    print
}' 1000G_Phase3_GRCh37.tsv \
> 1000G_Phase3_GRCh37_SNP.tsv