



#set directory to directory where bash script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR

bg_dir="/beevol/home/riemondy/Projects/rnastruct/data/mismatches/may2019/normalized_bg"

#for bg_file in $bg_dir/*.gz
#  do echo $bg_file
#      cp $bg_file .
#  done
 
chroms="/beevol/home/riemondy/Projects/rnastruct/pipeline/norm_enzymatic_data/chroms.txt"

for file in *.gz
do echo $file
  gunzip -c $file \
    | cut -f 1-4 \
    > tmp.bg

  bedGraphToBigWig \
      tmp.bg \
      $chroms \
      ${file/.bg.gz/.bb}
  rm tmp.bg
done

