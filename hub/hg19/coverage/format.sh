chroms="/beevol/home/riemondy/Projects/rnastruct/pipeline/norm_enzymatic_data/chroms.txt"

#for file in *pos*.bedgraph.gz
#do echo $file
#  gunzip -c $file \
#    | cut -f 1-4 \
#    > tmp.bg
#
#  bedGraphToBigWig \
#      tmp.bg \
#      $chroms \
#      ${file/.bedgraph.gz/.bb}
#  rm tmp.bg
#done

for file in *neg*.bedgraph.gz
do echo $file
  gunzip -c $file \
    | awk '{OFS=FS="\t"} {print $1,$2,$3,"-"$4}' \
    > tmp.bg

  bedGraphToBigWig \
      tmp.bg \
      $chroms \
      ${file/.bedgraph.gz/.bb}
  rm tmp.bg
done
