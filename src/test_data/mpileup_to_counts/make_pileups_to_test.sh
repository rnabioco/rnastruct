samtools mpileup \
  -f ../fasta/chr12_16_17_subset.fa.gz \
  -d 100000 \
  -L 100000 \
  -x -Q 0 -q 0 -B -A -C 0 \
  ../chr16_fus_dms.bam \
  | awk '($2 >= 31199140 && $2 <= 31199149) || ($2 >= 31199180 && $2 <=  31199189) || ($2 >= 31199196 && $2 <= 31199200)' \
      > pileups_to_test.txt

