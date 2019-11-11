#! /usr/bin/env bash 

c_flag='false'
l_flag='false'
verbose='false'

print_usage() {
  printf "Usage: run_tests.sh [-c 'clean up test data files']\n"
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

set -x

if [[ "$l_flag" == "true" ]]; then
    echo "running local tests"
    
  fa=$HOME/Projects/shared_dbases/genomes/human/chr_appended_grch37/Homo_sapiens.GRCh37.dna.primary_assembly_chrappended.fa 
  
  base_counter.py \
      -b test_data/small.fwd.bam \
      -f $fa \
      -n 'ACTG' \
      -t 2 \
      -o "test_data/base_counter/fwd_" \
      -L 'paired' \
      -d 100 \
      -v
  
  base_counter.py \
      -b test_data/small.rev.bam \
      -f $fa \
      -n 'ACTG' \
      -t 2 \
      -o "test_data/base_counter/rev_" \
      -L 'paired' \
      -d 100 \
      -v
  
  base_counter.py \
      -b test_data/small.bam \
      -f $fa \
      -n 'ACTG' \
      -t 2 \
      -o "test_data/base_counter/both_" \
      -L 'paired' \
      -d 100 \
      -v
fi

pytest .

