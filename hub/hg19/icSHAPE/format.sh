
files=(
hek_ch_vitro_minus.ct
hek_ch_vitro_plus.ct
hek_ch_vivo_minus.ct
hek_ch_vivo_plus.ct
hek_cy_vitro_minus.ct
hek_cy_vitro_plus.ct
hek_cy_vivo_minus.ct
hek_cy_vivo_plus.ct
hek_np_vitro_minus.ct
hek_np_vitro_plus.ct
hek_np_vivo_minus.ct
hek_np_vivo_plus.ct
)

chain="../rnastructurome/hg38ToHg19.over.chain.gz"
chroms="../rnastructurome/hg19.chrom.sizes"

for file in ${files[@]}
do echo $file
  awk 'NR > 1' $file | \
    liftOver stdin $chain stdout ${file/.ct/.unmapped} | \
    sort -k1,1 -k2,2n -k3,3n | \
    bedtools merge -i - -d -1 -c 4 -o mean \
    > ${file/.ct/.hg19}
  
  bedGraphToBigWig \
    ${file/.ct/.hg19} \
    $chroms \
    ${file/.ct/.bw}

  rm ${file/.ct/.hg19}

done
