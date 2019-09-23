#! /usr/bin/env bash
#BSUB -n 1
#BSUB -J rnafold
#BSUB -e err_%J.txt
#BSUB -o out_%J.txt
#BSUB -R "rusage[mem=4] select[mem>4] span[hosts=1]"
#BSUB -q rna

. /usr/share/Modules/init/bash
module load modules modules-init modules-python
module load viennarna

cd no_slop
outdir="RNAfold"
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
#  files=($depth"_cfiles/"*.constraints)
#  
#  outfile=$outdir"/"$depth"_rnafold.fasta"
#  rm -f $outfile
#  
#  for file in "${files[@]}"
#  do echo 
#      RNAfold \
#      --noPS \
#      -C \
#      $file >> $outfile
#  done
#done
#
#
## no dms constraints
#depths=(
#10
#15
#50
#100
#)
#
#for depth in "${depths[@]}"
#do
#  files=($depth"_cfiles/"*.constraints)
#  
#  outfile=$outdir"/"$depth"_rnafold_nodms.fasta"
#  rm -f $outfile
#  
#  for file in "${files[@]}"
#  do echo 
#      RNAfold \
#      --noPS \
#      $file >> $outfile
#  done
#done
#
#
## use shape normalized values
#
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
#  outfile=$outdir"/"$depth"_rnafold_shape.fasta"
#  rm -f $outfile
#  
#  for file in "${files[@]}"
#  do echo 
#      python ~/Projects/rnastruct/bin/mapToFastaShape.py \
#        --uConvert \
#        $file
#
#      RNAfold \
#      --noPS \
#      --shape=${file/.map/.shape} \
#      --shapeMethod="D" \
#      ${file/.map/.fa} >> $outfile
#  done
#done



# use approach from Wang et al 2019
# (https://doi.org/10.1016/j.ymeth.2018.11.018) 

depths=(
10
15
50
100
)

for depth in "${depths[@]}"
do
  files=($depth"_mapfiles/"*.map)
  
  outfile=$outdir"/"$depth"_rnafold_wangetal.fasta"
  rm -f $outfile
  
  for file in "${files[@]}"
  do echo 
      RNApvmin \
          -d2 \
          --shapeConversion=C0.03 \
          --tauSigmaRatio=1.0 \
          --minimizer=default \
          --initialVector=0 \
          ${file/.map/.shape} \
          < ${file/.map/.fa} \
          > ${file/.map/.sc} 

      RNAfold \
      --noPS \
      --shape=${file/.map/.sc} \
      --shapeMethod="W" \
      ${file/.map/.fa} >> $outfile
  done
done

