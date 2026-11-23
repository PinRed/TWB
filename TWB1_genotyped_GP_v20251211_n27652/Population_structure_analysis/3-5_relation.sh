ARRAY=TWB1
EAS=${ARRAY}_EAS.id
relation=${ARRAY}_filter.relatedsamples

awk 'NR==FNR {a[$1]=1; next} !($1 in a) {print $1, $1}' "$relation" "$EAS" > ${ARRAY}_unrelated_keep.id

rm -rf ${ARRAY}_S04_prune.prune.{in,out} \
    ${ARRAY}_S04_ldpruned.{bed,bim,fam} \
    ${ARRAY}_S04_qc.{bed,bim,fam} \
    ${ARRAY}_S04_EAS.{bed,bim,fam} \
    ${ARRAY}_relatedness_deg2plus.{kin0,kin0.zst,king} \
    ${ARRAY}_S04_duplicate.con