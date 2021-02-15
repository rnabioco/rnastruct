# rnastruct
[![Build Status](https://travis-ci.org/rnabioco/rnastruct.svg?branch=master)](https://travis-ci.org/rnabioco/rnastruct)

# Getting started

## Requirements
This pipeline requires the following programs and has been tested on
MacOS and linux:

[`samtools`](http://www.htslib.org/download/)  
[`python3`](https://www.python.org/downloads/)   

Also python packages `pandas`, `numpy`, `pysam`, and `dask`. See
`src/requirements.txt` for additional requirements

## Install 

```bash
git clone git@github.com:rnabioco/rnastruct
cd rnastruct/src
```

Add `src` directory to your path by editing your `.bash_profile` or `.bashrc`
```bash
export PATH="/Your/path/to/rnastruct/src:$PATH"
```


Test your install

```bash
$ base_counter.py -h
usage: base_counter.py [-h] -b BAM -f FASTA [-L LIBRARY] [-s STRANDEDNESS]
                       [-d DEPTH] [-l DELETION_LENGTH] [-o OUTPRE]
                       [-t THREADS] [-n NUCLEOTIDES] [-ss]
                       [-c CHROMS_TO_SKIP [CHROMS_TO_SKIP ...]] [-r REGION]
                       [-v] [-k] [--pileup-args KEY=VALUE [KEY=VALUE ...]]
                       [--pileup-arg-fn PILEUP_ARG_FN]

    Parse bam file and enumerate mismatches, insertions, and deletions per
    nucleotide. Generates bedgraphs for mismatches, insertions and
    deletions, as well as a summary table with aggregated information

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     indexed bam file input,
                                                by default the bam will be split into
                                                forward and reverse bams,
                                                pass the forward and reverse bams
                                                to bypass splitting the bams
                                                (i.e. forward.bam,reverse.bam)

  -f FASTA, --fasta FASTA
                        path to fasta file indexed with samtools faidx,
                        passed to samtools mpileup

  -L LIBRARY, --library LIBRARY
                        library type: either 'paired' or 'single' (default: paired)

  -s STRANDEDNESS, --strandedness STRANDEDNESS
                        strandedness of library:

                          'fr-firststrand' = first strand sequenced as R1 (default)
                             i.e. R1 alignments are the reverse complement of RNA,
                             for paired-end R2 alignments are the same sequence as the RNA

                          'fr-secondstrand' = second strand sequenced as R1
                             i.e. R1 alignments are the same sequence as the RNA
                             for paired-end R2 alignments are the reverse complement of the RNA

                          'unstranded' = report strandedness without respect to R1 or R2
                          (default: fr-firststrand)

  -d DEPTH, --depth DEPTH
                        minimum read coverage required for
                        reporting mismatch or indel frequencies (default: 5)

  -l DELETION_LENGTH, --deletion_length DELETION_LENGTH
                        do not count deletions greater than or equal
                        to deletion_length (default: 4)

  -o OUTPRE, --outpre OUTPRE
                        prefix for output files

  -t THREADS, --threads THREADS
                        Threads to use when running mpileup. If
                        threads is > 1 then the mpileup command will be
                        split up by chromosome to run using multiple
                        threads (default: 1)

  -n NUCLEOTIDES, --nucleotides NUCLEOTIDES
                        Nucleotides to use for computing
                        mismatch and indel ratios. Provide
                        as a string. i.e. to report
                        for all nucleotides. "ATCG"
                        (default: AC)

  -ss, --skip-single-reads
                        If set, skip single end reads when encountered
                        in a paired-end library.
                        (default: False , single end reads are
                        counted by default)

  -c CHROMS_TO_SKIP [CHROMS_TO_SKIP ...], --chroms-to-skip CHROMS_TO_SKIP [CHROMS_TO_SKIP ...]
                        space separated list of chroms to ignore
                        (default: None)

  -r REGION, --region REGION
                        region to query (as samtools region string)
                        (default: None)

  -v, --verbose         print run information (default: False)
  -k, --keep-temp-files
                        don't delete temp files (default: True)
  --pileup-args KEY=VALUE [KEY=VALUE ...]
                        Add key/value pileup arguments to overwrite
                        defaults, also can use --pileup-arg-fn to specify args
  --pileup-arg-fn PILEUP_ARG_FN
                        specify custom pileup.yaml file to overwrite default
                        arguments
```

## Generate mismatch and indel bedgraphs

`base_counter.py` is a python script that generates per nucleotide counts of indels and mismatches in 
a tabix indexed flat tsv file. 

The script generates and parses the output from running samtools
[`mpileup`](http://www.htslib.org/doc/samtools.html) and 
`mpileup_to_counts.py`

```bash
base_counter.py \
  -b alignments.bam \
  -f genome.fasta \
  -n 'ATCG' \
  -o "outfile/out_" 

```

The following files will be generated:

```
pileup_table.tsv.bgz
pileup_table.tsv.bgz.tbi
```

`pileup_table.tsv.bgz` contains per nucleotide counts of read depth, counts
of each nucleotide, indels, and deletions.


## Subtract background signals and normalize

```
filter_and_norm.py \
  -t dms_pileup_table.tsv.bgz \
  -u untreated_pileup_table.tsv.bgz \
  -p norm_pileup_table.tsv.bgz

```


## Generate bedgraphs from tabix pileup tables

```
tabix_to_bedgraph.py \
  -i norm_pileup_table.tsv.bgz \
  -c 20 # column to retrieve as score
  -s "+" #strand to retrieve

```

