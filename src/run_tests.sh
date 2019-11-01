#! /usr/bin/env bash 

fa=$HOME/Projects/shared_dbases/genomes/human/chr_appended_grch37/Homo_sapiens.GRCh37.dna.primary_assembly_chrappended.fa 

base_counter.py \
    -b test_data/small.fwd.bam \
    -f $fa \
    -n 'ACTG' \
    -t 2 \
    -o "test_data/base_counter/fwd_" \
    -L 'paired' \
    -d 100 


base_counter.py \
    -b test_data/small.rev.bam \
    -f $fa \
    -n 'ACTG' \
    -t 2 \
    -o "test_data/base_counter/rev_" \
    -L 'paired' \
    -d 100 

base_counter.py \
    -b test_data/small.bam \
    -f $fa \
    -n 'ACTG' \
    -t 2 \
    -o "test_data/base_counter/both_" \
    -L 'paired' \
    -d 100 

pytest .

