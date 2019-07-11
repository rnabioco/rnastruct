grep -v "^#" Parker_et_al_supp_data_file_1.tab \
  |  cut -f 2,3,4,5,6,7,8,9,10 \
  > evofam.bed

 
liftOver \
    -bedPlus=3 \
    -tab \
     evofam.bed \
     hg18ToHg19.over.chain.gz \
     evofam_hg19.bed \
     out.bed

sort -k1,1 -k2,2n evofam_hg19.bed > evofam_hg19.sorted.bed

bedToBigBed \
    -type=bed6+3 \
    -as=evofam.as \
    -tab \
    evofam_hg19.sorted.bed \
    http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes \
    evofam_hg19.bb
