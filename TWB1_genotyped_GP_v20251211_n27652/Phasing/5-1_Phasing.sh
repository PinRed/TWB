ARRAY="TWB1"
INPUT_DIR="../SNP_QC/Harmonization/VCF_Files"
OUT_DIR="Phased_VCF"
mkdir -p ${OUT_DIR}

for chr in {1..22}; do
    echo "Processing chromosome ${chr}"

    java -Xmx16g -jar ../../beagle/beagle.08Feb22.fa4.jar \
        gt=${INPUT_DIR}/${ARRAY}_chr${chr}.vcf.gz \
        out=${OUT_DIR}/${ARRAY}_chr${chr}_phased \
        nthreads=8

done