



#set directory to directory where bash script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

#wt="~/Projects/rnastruct/data/enzymatic_data_may2019/norm_data/data/struct_score/WT/single_files/WT.tsv.gz"
#c4="~/Projects/rnastruct/data/enzymatic_data_may2019/norm_data/data/struct_score/C4/single_files/C4.tsv.gz"

#cp $wt .
#cp $c4 .

chroms="/beevol/home/riemondy/Projects/rnastruct/pipeline/norm_enzymatic_data/chroms.txt"

for file in C4.tsv.gz
do echo $file
  table_to_bedgraph.py \
      -i $file \
      -b 2 \
      -e 3 \
      -sc 11 \
      -c 10 \
      -s "+" \
      -n \
      | sort -k1,1 -k2,2n  > $file".tmp_pos.bg"

  bedGraphToBigWig \
      $file".tmp_pos.bg" \
      $chroms \
      ${file/.tsv.gz/.pos.bw}
  rm $file".tmp_pos.bg"
  
  table_to_bedgraph.py \
      -i $file \
      -b 2 \
      -e 3 \
      -sc 11 \
      -c 10 \
      -s "-" \
      -n \
      | sort -k1,1 -k2,2n  > $file".tmp_neg.bg"

  bedGraphToBigWig \
      $file".tmp_neg.bg" \
      $chroms \
      ${file/.tsv.gz/.neg.bw}
  rm $file".tmp_neg.bg"
done

