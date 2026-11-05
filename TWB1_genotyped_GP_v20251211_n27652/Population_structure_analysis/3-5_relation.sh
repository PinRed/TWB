ARRAY=TWB1
EAS=${ARRAY}_EAS.id
relation=${ARRAY}_filter.relatedsamples

awk 'NR==FNR {a[$1]=1; next} !($1 in a) {print $1, $1}' "$relation" "$EAS" > ${ARRAY}_unrelated_keep.id