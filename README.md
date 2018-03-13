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

Also `samtools` needs to be installed and in your path

Test your install
```
$ count_mutations.py -h
usage: count_mutations.py [-h] -b BAM -f FASTA [-L LIBRARY] [-s STRANDEDNESS]
                          [-p PILEUP] [-d DEPTH] [-l DELETION_LENGTH]
                          [-o OUTPRE] [-t THREADS] [-v VERBOSE]

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
                        library type: either 'paired' or 'single' (default: single)

  -s STRANDEDNESS, --strandedness STRANDEDNESS
                        strandedness of library:

                          'fr-firststrand' = first strand sequenced as R1 (default)
                             i.e. R1 alignments are the reverse complement of RNA,
                             for paired-end R2 alignments are the same sequence as the RNA

                          'fr-secondstrand' = second strand sequenced as R1
                             i.e. R1 alignmnets are the same sequence as the RNA
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
                        do not count orphan reads (paired end)
                        do not double count if paired end reads overlap

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

  -v VERBOSE, --verbose VERBOSE
                        print run information (default: False)
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

The following files will be generated:

```
pileup_table.tsv.gz
neg_depth.bedgraph.gz
neg_deletions.bedgraph.gz
neg_insertions.bedgraph.gz
neg_mismatches.bedgraph.gz
pos_depth.bedgraph.gz
pos_deletions.bedgraph.gz
pos_insertions.bedgraph.gz
pos_mismatches.bedgraph.gz
```

`pileup_table.tsv.gz` contains per nucleotide counts of read depth, counts
of each nucleotide, insertion, and deletion per strand.

The `.bedgraph.gz` files are strand specific bedgraphs for mismatches and
indels. 


