# Structure score 

Structure scores are calculated using the output from running `MultiBamSummary`
from the `Deeptools` package on the bam files from ss-RNase or ds-RNase
treated samples.  


e.g.
```bash
multiBamSummary \
  bins \
  -b dsrnase.bam ssrnase.bam \
  -bs 1 \
  -v \
  -out output.npz 
  --outRawCounts output.tsv

# remove regions with no counts and sort
awk 'NR > 1 && ($4 + $5) != 0' output.tsv \
  | sort -k1,1 -k2,2n -k3,3n - \
  | gzip > output.tsv.gz
```


```bash
python structure_score.py \
  -c output.tsv.gv \ # text file from Deeptools
  -m "counts" \ # count based normalization (or coverage)
  -r 125456,432135 # library sizes 
```

The output column are 
chrom, start (0-based), end, dsrnase_counts, srnase_counts, norm_dsrnase,
norm_ssrnase, structure_score

