WORKING_PREFIX=../GCR/TWB1_S01_GCR
ARRAY=TWB1

plink2 \
  --bfile "${WORKING_PREFIX}" \
  --check-sex max-female-xf=0.2 min-male-xf=0.8 max-female-yrate=0.1 min-male-yrate=0.9 \
  --out ${ARRAY}_S02_sex

# View Mean
awk '
NR>1 {
    key=$4
    n[key]++
    if ($6!="NA") {
        sumF[key]+=$6
        nF[key]++
    }
    if ($7!="NA") {
        sumY[key]+=$7
        nY[key]++
    }
}
END {
    print "SNPSEX\tN\tMean_F\tMean_YRATE"
    for (x in n) {
        meanF=(nF[x] ? sumF[x]/nF[x] : "NA")
        meanY=(nY[x] ? sumY[x]/nY[x] : "NA")
        print x, n[x], meanF, meanY
    }
}' OFS='\t' ${ARRAY}_S02_sex.sexcheck | column -t -s $'\t'

# View table
awk '
BEGIN {
    FS=OFS="\t"
}
NR>1 {
    total++
    if ($5=="OK") ok++
    if ($5=="PROBLEM") problem++
    if ($3=="0") missing_recorded++
    if ($3=="1" && $4=="2") documented_mismatch++
    if ($3=="2" && $4=="1") documented_mismatch++
}
END {
    print "Metric\tCount"
    print "Total_samples\t" total
    print "STATUS_OK\t" ok
    print "STATUS_PROBLEM\t" problem
    print "Missing_recorded_sex\t" missing_recorded
    print "Confirmed_documented_vs_genetic_mismatch\t" documented_mismatch
}' ${ARRAY}_S02_sex.sexcheck > ${ARRAY}_S02_sex.summary.tsv

awk 'BEGIN {FS=OFS="\t"} NR>1 && $5=="PROBLEM" {print $1, $2}' ${ARRAY}_S02_sex.sexcheck > ${ARRAY}_S02_sex.problems.txt