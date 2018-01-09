# rnastruct

# Getting started

## Requirements
This pipeline requires the following programs and has been tested on
MacOS:

[`samtools`](http://www.htslib.org/download/)  
[`python3`](https://www.python.org/downloads/)   
[`pandas`](https://pandas.pydata.org/) python package  
A C/C++ compiler with C++11 support (i.e. gcc or clang, available for
macOS by installing the `Xcode` command line tools)

## Install 

```bash
git clone --recursive git@github.com:rnabioco/rnastruct
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
```
count_mutations.py -h
usage: count_mutations.py [-h] -b BAM -f FASTA [-p PILEUP] [-d DEPTH]
                          [-o OUTPRE] [-t THREADS] [-v VERBOSE]

Parse bam file and enumerate mismatches, insertions, and deletions per
nucleotide. Generates bedgraphs for mismatches, insertions and deletions, as
well as a summary table with aggregated information

optional arguments:
  -h, --help            show this help message and exit
  -b BAM, --bam BAM     indexed bam file input
  -f FASTA, --fasta FASTA
                        path to fasta file indexed with samtools faidx, passed
                        to samtools mpileup
  -p PILEUP, --pileup PILEUP
                        additional command line arguments to pass to samtools
                        mpileup, by default -f is set by the --fasta argument
                        to this script The following arguments are set by
                        default, but can be modified. --ff
                        UNMAP,SECONDARY,QCFAIL,DUP (filter alignments) -B
                        (disable BAQ calculation) -d 10000000 (use up to 1e7
                        reads per base) -L 10000000 (use up to 1e7 reads per
                        base for indel calling) do not count orphan reads
                        (paired end) do not double count if paired end reads
                        overlap
  -d DEPTH, --depth DEPTH
                        minimum read coverage required for reporting mismatch
                        or indel frequencies. Default = 5
  -o OUTPRE, --outpre OUTPRE
                        prefix for output files
  -t THREADS, --threads THREADS
                        Threads to use when running mpileup. If threads is > 1
                        then the mpileup command will be split up by
                        chromosome to run using multiple threads default = 1
  -v VERBOSE, --verbose VERBOSE
                        print run information
```

## Generate mismatch and indel bedgraphs

`count_mutations.py` is a python script that generates per nucleotide counts of indels and mismatches in bedgraph
format. 

The script generates and parses the output from running samtools
[`mpileup`](http://www.htslib.org/doc/samtools.html) and 
[`mpileup2ReadCounts`](https://github.com/kriemo/mpileup2readcounts).

```bash
count_mutations.py \
  -b alignments.bam \
  -f genome.fasta 
```

7 files will be generated:

```
pileup_table.tsv.gz
neg_deletions.bedgraph.gz
neg_insertions.bedgraph.gz
neg_mismatches.bedgraph.gz
pos_deletions.bedgraph.gz
pos_insertions.bedgraph.gz
pos_mismatches.bedgraph.gz
```

`pileup_table.tsv.gz` contains per nucleotide counts of read depth, counts
of each nucleotide, insertion, and deletion per strand.

The `.bedgraph.gz` files are strand specific bedgraphs for mismatches and
indels. 


