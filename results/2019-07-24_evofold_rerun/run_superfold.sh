#! /usr/bin/env bash
#BSUB -n 1
#BSUB -J convert[1-500]%100
#BSUB -e err_%J.txt
#BSUB -o out_%J.txt
#BSUB -R "rusage[mem=4] select[mem>4] span[hosts=1]"
#BSUB -q rna

. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load gcc

bin="/beevol/home/riemondy/Projects/rnastruct/bin/Superfold_1.0"

FILES=(*.map)

file=${FILES[$(($LSB_JOBINDEX - 1))]}

outdp=$(python $bin/Superfold.py \
    --np 1 \
    --noPVclient \
    $file) 

for ct in $outdp/*.ct
  do echo $ct
      ct2dot $ct 1 ${ct/.ct/.db}
  done
