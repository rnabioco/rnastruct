#! /usr/bin/env bash
#BSUB -n 1
#BSUB -J rnastructure 
#BSUB -e err_%J.txt
#BSUB -o out_%J.txt
#BSUB -R "rusage[mem=4] select[mem>4] span[hosts=1]"
#BSUB -q rna

. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load gcc

cd no_slop

outdir="RNAstructure"
mkdir -p $outdir

#depths=(
#10
#15
#50
#100
#)
#
#for depth in "${depths[@]}"
#do
#  files=($depth"_mapfiles/"*.map)
#  
#  outfile=$outdir"/"$depth"_rnastructure.db"
#  rm -f $outfile
#  
#  for file in "${files[@]}"
#  do echo $file 
#    
#    #python ~/Projects/rnastruct/bin/mapToFastaShape.py \
#    #    --uConvert \
#    #    $file
#    
#    Fold \
#        ${file/.map/.fa} \
#        - \
#        -k \
#        -dms ${file/.map/.shape} \
#        -md 600 \
#        -m 100 \
#        >> $outfile
#
#  done
#done


depths=(
10
15
50
100
)

for depth in "${depths[@]}"
do
  files=($depth"_mapfiles/no_dms/"*.map)
  
  outfile=$outdir"/"$depth"_rnastructure_nodms.db"
  rm -f $outfile

  for file in "${files[@]}"
  do echo $file 
    
    python ~/Projects/rnastruct/bin/mapToFastaShape.py \
        --uConvert \
        $file
    
    Fold \
        ${file/.map/.fa} \
        - \
        -k \
        -md 600 \
        -m 100 \
        >> $outfile
  done
done
