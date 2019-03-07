#! /usr/bin/env bash
#BSUB -J rnastruct 
#BSUB -o rnastruct_%J.out
#BSUB -e rnastruct_%J.err
#BSUB -R "select[mem>120] rusage[mem=120] span[hosts=1]" 
#BSUB -q rna
#BSUB -n 24

set -o nounset -o pipefail -o errexit -x

#### load necessary modules ####

. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load gcc
module load samtools
module load python/3.6.3

BAM="$HOME/Projects/rnastruct/data/newdata/WT_invivo_DMS_all_merge_sorted.bam"
FA="$HOME/Projects/shared_dbases/genomes/human/chr_appended_grch37/Homo_sapiens.GRCh37.dna.primary_assembly_chrappended.fa"
OUTDIR="$HOME/Projects/rnastruct/data/mismatches/new_data/"

./count_mutations.py \
  -f $FA \
  -b $BAM \
  -o $OUTDIR \
  -t 24 \
  -d 10 \
  -L 'paired' \
  -v 
