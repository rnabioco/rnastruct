chroms="/beevol/home/riemondy/Projects/rnastruct/pipeline/norm_enzymatic_data/chroms.txt"

for file in *neg*.bw
do echo $file

  bigWigToBedGraph $file stdout \
      | awk '{OFS=FS="\t"} {print $1,$2,$3,"-"$4}' \
    > tmp.bg

  bedGraphToBigWig \
      tmp.bg \
      $chroms \
      $file
  rm tmp.bg
done
