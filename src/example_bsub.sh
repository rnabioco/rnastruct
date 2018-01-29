#! /usr/bin/env bash
#BSUB -J rnastruct 
#BSUB -o rnastruct_%J.out
#BSUB -e rnastruct_%J.err
#BSUB -R "select[mem>10] rusage[mem=10] span[hosts=1]" 
#BSUB -q normal
#BSUB -n 12

set -o nounset -o pipefail -o errexit -x

#### load necessary modules ####

. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load gcc
module load samtools
module load python3

BAM="$HOME/Projects/rnastruct/tests/data/med.bam"
FA="$HOME/Projects/shared_dbases/genomes/GRCh38.primary_assembly.genome.fa"
OUTDIR="$HOME/Projects/rnastruct/data/mismatches/medium/"

count_mutations.py \
  -f $FA \
  -b $BAM \
  -o $OUTDIR
  