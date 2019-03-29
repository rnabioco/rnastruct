#! /usr/bin/env

bam=$1
outdir=$2

basebam=$(basename $bam)

mkdir -p $outdir

outname=$outdir"/"$basebam

# include reads that are 2nd in a pair (128);
# exclude reads that are mapped to the reverse strand (16)
samtools view -b -f 128 -F 16 $bam > $outname".fwd1.bam"

# exclude reads that are mapped to the reverse strand (16) and
# first in a pair (64): 64 + 16 = 80
samtools view -b -f 80 $bam > $outname".fwd2.bam"

# combine the temporary files
samtools merge -f $outname".fwd.bam" $outname".fwd1.bam" $outname".fwd2.bam"

# index the filtered BAM file
samtools index $outname".fwd.bam"

# include reads that map to the reverse strand (128)
# and are second in a pair (16): 128 + 16 = 144
samtools view -b -f 144 $bam > $outname".rev1.bam" 

# include reads that are first in a pair (64), but
# exclude those ones that map to the reverse strand (16)
samtools view -b -f 64 -F 16 $bam > $outname".rev2.bam"

# merge the temporary files
samtools merge -f $outname".rev.bam" $outname".rev1.bam" $outname".rev2.bam"

# index the merged, filtered BAM file
samtools index $outname".rev.bam"

rm $outname".fwd1.bam" $outname".fwd2.bam"
rm $outname".rev1.bam" $outname".rev2.bam"
