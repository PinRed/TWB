mkdir -p harmonized harmonized_common

plink2 \
    --bfile chr_pos/1kgp_pos \
    --extract chr_pos/common.pos \
    --make-bed \
    --out harmonized/1kgp

plink2 \
    --bfile chr_pos/sgdp_pos \
    --extract chr_pos/common.pos \
    --make-bed \
    --out harmonized/sgdp

plink2 \
    --bfile chr_pos/twb2_pos \
    --extract chr_pos/common.pos \
    --make-bed \
    --out harmonized/twb2


############################################################
# 4. 用 1KGP 當 allele 標準
#
#    判斷 SGDP / TWB2 是否：
#
#     已經相同
#     complement
#
#    只用 awk 產生 flip list
############################################################


############################################################
# 5. 正確產生 flip list
############################################################

awk '
function comp(a) {
    if (a=="A") return "T"
    if (a=="T") return "A"
    if (a=="C") return "G"
    if (a=="G") return "C"
    return a
}

function is_ambiguous(a1, a2) {
    return ((a1=="A" && a2=="T") || (a1=="T" && a2=="A") || (a1=="C" && a2=="G") || (a1=="G" && a2=="C"))
}

FNR==NR {
    ref_a1[$2] = $5
    ref_a2[$2] = $6
    next
}

{
    id=$2
    a1=$5
    a2=$6

    r1=ref_a1[id]
    r2=ref_a2[id]

    # 過濾 non Strand Ambiguous SNPs
    if (is_ambiguous(r1, r2) || is_ambiguous(a1, a2)) {
        print id > mismatch
        next
    }

    # 完全相同
    if (a1==r1 && a2==r2)
        next

    # complement
    if (comp(a1)==r1 && comp(a2)==r2) {
        print id
        next
    }


    # 其他就是 allele mismatch
    print id > mismatch
}

' mismatch="harmonized/sgdp.mismatch" \
  harmonized/1kgp.bim \
  harmonized/sgdp.bim \
  > harmonized/sgdp.flip


awk '
function comp(a) {
    if (a=="A") return "T"
    if (a=="T") return "A"
    if (a=="C") return "G"
    if (a=="G") return "C"
    return a
}

function is_ambiguous(a1, a2) {
    return ((a1=="A" && a2=="T") || (a1=="T" && a2=="A") || (a1=="C" && a2=="G") || (a1=="G" && a2=="C"))
}

FNR==NR {
    ref_a1[$2] = $5
    ref_a2[$2] = $6
    next
}

{
    id=$2
    a1=$5
    a2=$6

    r1=ref_a1[id]
    r2=ref_a2[id]

    # 過濾 non Strand Ambiguous SNPs
    if (is_ambiguous(r1, r2) || is_ambiguous(a1, a2)) {
        print id > mismatch
        next
    }

    if (a1==r1 && a2==r2)
        next


    if (comp(a1)==r1 && comp(a2)==r2) {
        print id
        next
    }

    print id > mismatch
}

' mismatch="harmonized/twb2.mismatch" \
  harmonized/1kgp.bim \
  harmonized/twb2.bim \
  > harmonized/twb2.flip


############################################################
# 6. 把 SGDP / TWB2 的 strand flip
############################################################

plink2 \
    --bfile harmonized/sgdp \
    --flip harmonized/sgdp.flip \
    --make-bed \
    --out harmonized/sgdp_flip

plink2 \
    --bfile harmonized/twb2 \
    --flip harmonized/twb2.flip \
    --make-bed \
    --out harmonized/twb2_flip


############################################################
# 7. 只保留三個 dataset allele 都相容的 SNP
############################################################

awk '
BEGIN { OFS="\t" }

FNR==NR {
    ref1[$2]=$5
    ref2[$2]=$6
    next
}

{
    id=$2
    a1=$5
    a2=$6

    r1=ref1[id]
    r2=ref2[id]

    if (a1==r1 && a2==r2) 
        print id
}
' harmonized/1kgp.bim harmonized/sgdp_flip.bim \
    > harmonized/sgdp.compatible


awk '
BEGIN { OFS="\t" }

FNR==NR {
    ref1[$2]=$5
    ref2[$2]=$6
    next
}

{
    id=$2
    a1=$5
    a2=$6

    r1=ref1[id]
    r2=ref2[id]

    if (a1==r1 && a2==r2)
        print id
}
' harmonized/1kgp.bim harmonized/twb2_flip.bim \
    > harmonized/twb2.compatible


############################################################
# 8. 三者取交集
############################################################

sort -u harmonized/sgdp.compatible > harmonized/sgdp.compatible.sorted
sort -u harmonized/twb2.compatible > harmonized/twb2.compatible.sorted

comm -12 harmonized/sgdp.compatible.sorted \
         harmonized/twb2.compatible.sorted \
         > harmonized/final_common.pos


############################################################
# 9. 最後全部限制到 final_common.pos
############################################################

plink2 \
    --bfile harmonized/1kgp \
    --extract harmonized/final_common.pos \
    --make-bed \
    --out harmonized/1kgp_finalset

plink2 \
    --bfile harmonized/sgdp_flip \
    --extract harmonized/final_common.pos \
    --make-bed \
    --out harmonized/sgdp_finalset

plink2 \
    --bfile harmonized/twb2_flip \
    --extract harmonized/final_common.pos \
    --make-bed \
    --out harmonized/twb2_finalset



############################################################
# 10. 最後 ID = CHR:POS:A1:A2
############################################################

plink2 \
    --bfile harmonized/1kgp_finalset \
    --set-all-var-ids '@:#:$r:$a' \
    --make-bed \
    --out harmonized_common/1kgp_common

plink2 \
    --bfile harmonized/sgdp_finalset \
    --set-all-var-ids '@:#:$r:$a' \
    --make-bed \
    --out harmonized_common/sgdp_common

plink2 \
    --bfile harmonized/twb2_finalset \
    --set-all-var-ids '@:#:$r:$a' \
    --make-bed \
    --out harmonized_common/twb2_common


############################################################
# 11. 檢查三個 BIM 是否完全一致
############################################################

echo "=== SNP counts ==="

wc -l \
    harmonized_common/1kgp_common.bim \
    harmonized_common/sgdp_common.bim \
    harmonized_common/twb2_common.bim


echo "=== ID differences ==="

cut -f2 harmonized_common/1kgp_common.bim | sort > harmonized/1kgp.ids
cut -f2 harmonized_common/sgdp_common.bim | sort > harmonized/sgdp.ids
cut -f2 harmonized_common/twb2_common.bim | sort > harmonized/twb2.ids

diff harmonized/1kgp.ids harmonized/sgdp.ids || true
diff harmonized/1kgp.ids harmonized/twb2.ids || true


echo "=== Allele differences ==="

paste \
    harmonized_common/1kgp_common.bim \
    harmonized_common/sgdp_common.bim \
    harmonized_common/twb2_common.bim |
awk '
{
    if ($5 != $11 || $6 != $12 ||
        $5 != $17 || $6 != $18)
        print
}' | head -20

rm -rf harmonized 
