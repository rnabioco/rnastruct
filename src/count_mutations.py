#!/usr/bin/env python3
import argparse
import os 
import sys
import subprocess
import gzip
import tempfile
from collections import Counter


def generate_pileup(bam, fasta, additional_args):
    """ returns tempfile pointing to pileup output 
    
    args:
        bam = path to bam
        fasta = path to fasta
        additional_args = list of additional arguments for samtools pileup
    """


    temp_output = tempfile.NamedTemporaryFile(delete=False)

    pileup_cmd = "samtools " + \
                  "mpileup " + \
                  "-f " + fasta + " " + \
                  additional_args + " " + \
                  bam + " " + \
                  " | " + \
                  "./mpileup2readcounts/bin/mpileup2readcounts - "
                 
    print("pileup command is:\n" + pileup_cmd, file = sys.stderr)
    
    pileup_run = subprocess.run(pileup_cmd, 
            shell=True, stderr = sys.stderr, stdout = temp_output)
    
    return temp_output

def format_bedgraphs(fname):

    f = open(fname)
    mismatches = open("mismatch.bg", "w") 
    inserts = open("insertions.bg", "w")
    dels = open("deletions.bg", "w")
   
    header = f.readline()
    for line in f:
        line = line.rstrip()
        chrom, end, strand, depth, ref_base,\
        ref_count, acount, ccount, tcount, gcount, \
        ncount, delcount, inscount = line.split("\t")
        
        # one-based input, zero-based output
        start = int(end) - 1
        shared_cols = [chrom, start, end]
        
        mm_counts = [acount, ccount, tcount, gcount]
        mm_counts = [int(x) for x in mm_counts]

        if any([x > 0 for x in mm_counts]):
            mismatch = shared_cols + mm_counts 
            mismatch = [str(x) for x in mismatch]
            mismatches.write("\t".join(mismatch) + "\n")
        
        if int(inscount) > 0:
            insertions = shared_cols + [inscount]
            insertions = [str(x) for x in insertions]
            inserts.write("\t".join(insertions) + "\n")
        
        if int(delcount) > 0:
            deletions = shared_cols + [delcount]
            deletions = [str(x) for x in deletions]
            dels.write("\t".join(deletions) + "\n")
    

def generate_mismatch_profile(bam, fasta, additional_args):
    
    # generate per nucleotide mismatch and indel counts
    temp_output = generate_pileup(bam, fasta, additional_args)
    temp_output.close()
    
    # parse into bedgraphs
    format_bedgraphs(temp_output.name) 
    os.unlink(temp_output.name)

def main():
    
    parser = argparse.ArgumentParser(description="""
    Parse bam file and enumerate mismatches, insertions, and deletions per
    nucleotide. Generates bedgraphs for mismatches, insertions and
    deletions, as well as a summary table with aggregated information""")

    parser.add_argument('-b',
                        '--bam',
                        help ='indexed bam file input',
                        required = True)
    parser.add_argument('-f',
                        '--fasta',
                        help = """path to fasta file indexed with samtools
                        faidx, passed to samtools mpileup""",
                        required = True)
    parser.add_argument('-p',
                        '--pileup',
                        help = """additional command line arguments
                        to pass to samtools mpileup, by default
                        -f is set by the --fasta argument to this script
                        
                        The following arguments are set by default, but can
                        be modified.
                        --ff UNMAP,SECONDARY,QCFAIL,DUP (filter alignments)
                        -B (disable BAQ calculation) 
                        -d 10000000 (use up to 1e7 reads per base)
                        -L 10000000 (use up to 1e7 reads per base for indel
                        calling)
                        do not count orphan reads (paired end)
                        do not double count if paired end reads overlap
                        """, 
                        required = False,
                        default = " -d 10000000 -L 10000000 -B ")

    args=parser.parse_args()
    
    bam_name = args.bam
    pileup_args = args.pileup    
    fasta_name = args.fasta
    generate_mismatch_profile(bam_name, fasta_name, pileup_args) 
                
if __name__ == '__main__': main()

