# rnastruct
[![Build Status](https://travis-ci.org/rnabioco/rnastruct.svg?branch=master)](https://travis-ci.org/rnabioco/rnastruct)

# Getting started

## Requirements
This pipeline requires the following programs and has been tested on
MacOS:

[`samtools`](http://www.htslib.org/download/)  
[`bcftools`](http://www.htslib.org/download/)  
[`python3`](https://www.python.org/downloads/)   

Also python packages `pandas`, `cyvcf2`, `pysam`, and `dask`
A C/C++ compiler with C++11 support (i.e. gcc or clang, available for
macOS by installing the `Xcode` command line tools)

## Install 

```bash
git clone git@github.com:rnabioco/rnastruct
cd rnastruct/src
make
```

Add `src` directory to your path by editing your `.bash_profile`
```bash
export PATH="/Your/path/to/rnastruct/src:$PATH"
```

Source your `.bash_profile` 

```bash
source ~/.bash_profile
```


Test your install

```bash
$ base_counter.py -h
usage: base_counter.py [-h] -b BAM -f FASTA [-L LIBRARY] [-s STRANDEDNESS]
                       [-p PILEUP] [-d DEPTH] [-l DELETION_LENGTH] [-o OUTPRE]
                       [-t THREADS] [-n NUCLEOTIDES] [-ss] [-i]
                       [-c CHROMS_TO_SKIP [CHROMS_TO_SKIP ...]] [-v] [-k]
                       [--write-bgs]
                       [--default-pileup-args DEFAULT_PILEUP_ARGS]

    Parse bam file and enumerate mismatches, insertions, and deletions per
    nucleotide. Generates bedgraphs for mismatches, insertions and
    deletions, as well as a summary table with aggregated information

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     indexed bam file input

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

  -p PILEUP, --pileup PILEUP
                        additional command line arguments to pass to samtools mpileup
                        by default -f is set by the --fasta argument to this script

                        The following arguments are set by default, but can be modified.
                        --ff UNMAP,SECONDARY,QCFAIL,DUP (filter alignments)
                        -B (disable BAQ calculation)
                        -d 1000000 (use up to 1e6 reads per base)
                        -L 1000000 (use up to 1e6 reads per base for indel calling)
                        -A count orphan reads (paired end)
                        -x disable read-pair overlap detection

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

  -i, --tabix-index     If set, report pileup table
                        as bgzip'd and tabix indexed
                        (default: False)

  -c CHROMS_TO_SKIP [CHROMS_TO_SKIP ...], --chroms-to-skip CHROMS_TO_SKIP [CHROMS_TO_SKIP ...]
                        space separated list of chroms to ignore
                        (default: None)

  -v, --verbose         print run information (default: False)
  -k, --keep-temp-files
                        don't delete temp files (default: True)
  --default-pileup-args DEFAULT_PILEUP_ARGS
                        The following arguments are set by default
                        --ff UNMAP,SECONDARY,QCFAIL,DUP (filter alignments)
                        -B (disable BAQ calculation)
                        -d 1000000 (use up to 1e6 reads per base)
                        -I dont call indels
                        -A count orphan reads (paired end)
                        -x disable read-pair overlap detection
                        -a AD
                        pass a string to replace these args
                        (default:  --ff UNMAP,SECONDARY,QCFAIL,DUP -a AD -A -x -d 100000 -L 100000 -B -O v )
```

## Generate mismatch and indel bedgraphs

`base_counter.py` is a python script that generates per nucleotide counts of indels and mismatches in 
a tabix indexed flat tsv file. 

The script generates and parses the output from running bcftools
[`mpileup`](http://www.htslib.org/doc/bcftools.html) and 
`pileup_to_counts.py`

```bash
base_counter.py \
  -b alignments.bam \
  -f genome.fasta 
```

The following files will be generated:

```
pileup_table.tsv.bgz
pileup_table.tsv.bgz.tbi

neg_depth.bedgraph.gz
neg_indels.bedgraph.gz
neg_mismatches.bedgraph.gz
neg_mutations.bedgraph.gz
pos_depth.bedgraph.gz
pos_indels.bedgraph.gz
pos_mismatches.bedgraph.gz
pos_mutations.bedgraph.gz
```

`pileup_table.tsv.bgz` contains per nucleotide counts of read depth, counts
of each nucleotide, indels, and deletions.

The `.bedgraph.gz` files are strand specific bedgraphs for mismatches and
indels. 


