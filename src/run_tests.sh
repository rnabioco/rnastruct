#! /usr/bin/env bash 

c_flag='false'
l_flag='false'
verbose='false'

print_usage() {
  printf "Usage: run_tests.sh\n"
  printf " [-c 'clean up test data files']\n"
  printf " [-l 'run local tests']\n"
}

while getopts 'clv' flag; do
  case "${flag}" in
    c) c_flag='true' ;;
    l) l_flag='true' ;;
    v) verbose='true' ;;
    *) print_usage
       exit 1 ;;
  esac
done

if [[ "$c_flag" == 'true' ]]; then
   echo "removing preprocessed files"
   rm -rf test_data/base_counter
fi

set -e -x

if [[ "$l_flag" == "true" ]]; then
    echo "running local tests"
    
  fa="test_data/fasta/chr12_16_17_subset.fa.gz"
  data_dir="test_data"
  out_dir=$data_dir"/base_counter"
  bams=(
  small.fwd.bam
  small.rev.bam
  small.bam
  chr16_fus_dms.bam
  chr16_fus_ut.bam
  )

  ids=(
  fwd_
  rev_
  both_
  fus_dms_
  fus_ut_
  )
  
  vals=(
  100
  100
  100
  1
  1
  )
  
  for (( i=0; i<${#bams[@]}; i++ ))
  do 
    base_counter.py \
        -b $data_dir"/"${bams[$i]} \
        -f $fa \
        -n 'ACTG' \
        -t 2 \
        -o $out_dir"/"${ids[$i]} \
        -L 'paired' \
        -d ${vals[$i]} \
        -v
  done

  # filter and normalize
  filter_and_norm.py \
    -t $out_dir"/fus_dms_pileup_table.tsv.bgz" \
    -u $out_dir"/fus_ut_pileup_table.tsv.bgz" \
    -o $out_dir"/fus_norm_" \
    -d 1 \
    -p 2 \
    -tr 1.0 \
    -ur 0.02 \
    -n 'ATCG'

  # generate bedgraphs
  to_get=(
  5
  11
  12
  15
  17
  19
  21
  )
  cols=(
  depth_
  mmcount_
  indel_count_
  mutation_ratio_
  mutation_ratio_ut_
  mutation_ratio_bgsub_
  mutation_ratio_norm_
  )

  out_bg_dir=$out_dir"/bedgraph"
  mkdir -p $out_bg_dir

  for (( i=0; i<${#to_get[@]}; i++ ))
  do
    tabix_to_bedgraph.py \
        -i $out_dir"/fus_norm_pileup_table.tsv.bgz"\
        -c ${to_get[$i]} \
        -s "+" \
        > $out_bg_dir"/"${cols[$i]}"pos.bedgraph"
  done
fi

pytest .

