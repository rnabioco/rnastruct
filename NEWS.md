
## 2019-02-24 

  - Try to run segway.
  - Use the following file types:
    - Normalized reactivity (WT_genome_score_sorted_fixed.bedgraph)
      - convert 4th base to 0
    - Depth (parse from WT_genome_merge_score.bed) 
    - A-to-I editing (Use WT samples, need to average replicates,)
      - just editing frequencies only
      - start with just pos strand at this point
      - exclude chrX and chrY from data to match depth and reactivity
    - make a chromsizes fix (without chrX and chrY)

## 2018-05-11

  - Compute sensitivity versus specificity for 18srRNA structure
  - Try 2-pass alignment via Hisat2 followed by bowtie2
  - Figure out which regions:
    a) Have enough coverage to predict structure
    b) Have differences in reactivity between fast, normal, and slow
    polymerases
  - Add summary statistics for alignments (% mismatches, % aligned, etc.)


