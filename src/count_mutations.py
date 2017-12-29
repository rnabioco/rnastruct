#!/usr/bin/env python3
import argparse
import os 
import sys
import subprocess
import gzip
import tempfile
import pandas as pd
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
                  "./mpileup2readcounts/bin/mpileup2readcounts - " + \
                  " | gzip "
                 
    print("pileup command is:\n" + pileup_cmd, file = sys.stderr)
    
    print("formatted pileup output is here:\n" + temp_output.name, 
            file = sys.stderr)

    pileup_run = subprocess.run(pileup_cmd, 
            shell=True, stderr = sys.stderr, stdout = temp_output)
    
    print("completed pileup", file = sys.stderr)

    return temp_output

def format_bedgraphs(df, depth, prefix):
    """ take pandas dataframe and generate bedgraphs """

    df = df.assign(mismatch_ratio = lambda df:1 - (df.refcount / df.depth))
    df = df.assign(insertion_ratio = lambda df: df.inscount / df.depth)
    df = df.assign(deletion_ratio = lambda df: df.delcount / df.depth)
    
    df = df.assign(start = lambda df:df.pos - 1)
    df = df.rename(columns = {'pos':'end'})
    
    df = df[df.depth >= depth]

    df_mm = df[['chr', 'start', 'end', 'mismatch_ratio']]
    df_ins = df[['chr', 'start', 'end', 'insertion_ratio']]
    df_del = df[['chr', 'start', 'end', 'deletion_ratio']]

    df_mm = df_mm.sort_values(['chr', 'start'], ascending=[True, True])
    df_ins = df_ins.sort_values(['chr', 'start'], ascending=[True, True])
    df_del = df_del.sort_values(['chr', 'start'], ascending=[True, True])
    
    df_mm.to_csv(prefix + "mismatches.bg.gz", sep= "\t", index=False, header=False, compression='gzip')
    df_ins.to_csv(prefix + "insertions.bg.gz", sep= "\t", index=False, header=False, compression='gzip')
    df_del.to_csv(prefix + "deletions.bg.gz", sep= "\t", index=False, header=False, compression='gzip')

    
def parse_library_type(bam, libtype):
    pass

    
##    # 1. alignments of the second in pair if they map to the forward strand
#    # 2. alignments of the first in pair if their mate maps to the forward strand
#    
#    samtools view -b -f 128 -F 16 $DATA > fwd1.bam
#    samtools index fwd1.bam
#
#    samtools view -b -f 64 -F 32 $DATA > fwd2.bam
#    samtools index fwd2.bam
#    
#    samtools merge -f fwd.bam fwd1.bam fwd2.bam
#    samtools index fwd.bam
#
#    # 1. alignments of the second in pair if it maps to the reverse strand
#    # 2. alignments of the first in pair if their mates map to the reverse strand
#    
#    samtools view -b -f 144 $DATA > rev1.bam
#    samtools index rev1.bam
#
#    samtools view -b -f 96 $DATA > rev2.bam
#    samtools index rev2.bam
#
#    samtools merge -f rev.bam rev1.bam rev2.bam
#    samtools index rev.bam

def generate_mismatch_profile(bam, fasta, additional_args, depth):
    
    # generate per nucleotide mismatch and indel counts
    temp_output = generate_pileup(bam, fasta, additional_args)
    temp_output.close()
    
    # parse into bedgraphs (pos and neg)
    df = pd.read_table(temp_output.name,
              compression='gzip')
    
    df_pos = df[df.strand == "+"]
    df_neg = df[df.strand == "-"]

    format_bedgraphs(df_pos, depth, "pos_") 
    format_bedgraphs(df_neg, depth, "neg_")
    
    # cleanup
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
    
    parser.add_argument('-d',
                        '--depth',
                        help = """minimum read coverage required for
                        reporting mismatch or indel frequencies. 
                        Default = 5""", 
                        required = False, 
                        default = 5, type = float)

    args=parser.parse_args()
    
    bam_name = args.bam
    pileup_args = args.pileup    
    fasta_name = args.fasta
    depth = args.depth

    generate_mismatch_profile(bam_name, 
            fasta_name, 
            pileup_args,
            depth) 
                
if __name__ == '__main__': main()

