#!/usr/bin/env python3
import argparse
import os 
import sys
import subprocess
import gzip
import pandas as pd
import uuid 
import shutil
from datetime import datetime
from collections import Counter
from multiprocessing import Pool
from functools import partial

def retrieve_header(bam):
    """
    args:
      bam = path to bam
    
    return:
      list of contigs
    """
    get_header = ["samtools", "view", "-H", bam]

    header_call = subprocess.run(get_header, 
            stdout = subprocess.PIPE, 
            universal_newlines = True).stdout
    contigs = []
    for line in header_call.splitlines():
       
        # keep contig lines
        if not line.startswith("@SQ"):
            continue
        
        contig_id = line.split("\t")[1]
        # drop "SN" 
        if contig_id.startswith('SN:'):
            contig_id = contig_id[3:]
        else:
            print("non-standard contig found in header {}".format(contig_id, 
                        file = sys.stderr))

        contigs.append(contig_id)
    
    return contigs

def generate_pileup(bam, fasta, min_depth, outpre, additional_args, verbose = False):
    """ returns fileobject to pileup output 
    
    args:
        bam = path to bam
        fasta = path to fasta
        additional_args = list of additional arguments for samtools pileup
    """


    output = open(outpre + "pileup_table.tsv.gz", "w")

    pileup_cmd = "samtools " + \
                  "mpileup " + \
                  "-f " + fasta + " " + \
                  additional_args + " " + \
                  bam + " " + \
                  " | " + \
                  "mpileupToReadCounts -d " + \
                  str(min_depth) + \
                  " | gzip "
                 
    if verbose:

        pileup_run = subprocess.run(pileup_cmd, 
            shell=True, stderr = sys.stderr, stdout = output)
        
        print("pileup command is:\n" + pileup_cmd, file = sys.stderr)
        
        print("formatted pileup output is here:\n" + output.name, 
                file = sys.stderr)
        print("completed pileup {}".format(str(datetime.now())), 
                file = sys.stderr)
    else:
        pileup_run = subprocess.run(pileup_cmd, 
            shell=True, stderr = subprocess.PIPE, stdout = output)

    output.close()
    return output.name

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
    
    df_mm.to_csv(prefix + "mismatches.bedgraph.gz", sep= "\t", index=False, header=False, compression='gzip')
    df_ins.to_csv(prefix + "insertions.bedgraph.gz", sep= "\t", index=False, header=False, compression='gzip')
    df_del.to_csv(prefix + "deletions.bedgraph.gz", sep= "\t", index=False, header=False, compression='gzip')

    
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

def generate_mismatch_profile(input_bam, fasta, additional_args, depth, outpre,
        threads = 1, debug = False):
    
    if debug:
        print("started processing {}".format(str(datetime.now())),
            file = sys.stderr)

    # generate per nucleotide mismatch and indel counts
    if threads == 1:
        output = generate_pileup(input_bam, fasta, depth, outpre, additional_args)
    
    else:
        """ if multiple threads then run mpileup on each chromosome using
        the region arg. Write output to temp folder, and combine results
        using pd.concat(). Use concat instead of more traditional
        approaches as there is a header line for each output that would
        need to be dropped before combining.
        """
        contigs = retrieve_header(input_bam)
        
        if "-r " in additional_args:
            print("-r option is not allowed when running with multiple threads", file = sys.stderr)
        
        # generate list of new regional arguments to pass in parallel
        new_args = [additional_args + " -r " + x + " " for x in contigs]

        pool = Pool(threads)
        
        # build function obj
        func = partial(generate_pileup,
            input_bam,
            fasta,
            depth,
            verbose = debug)
        
        # make tmp directory
        tmp_dir = "tmpfiles-" + str(uuid.uuid4())
        if not os.path.exists(tmp_dir): os.makedirs(tmp_dir)
        
        # generate prefixes that with tmp_dir and contig id
        new_pres = [os.path.join(tmp_dir, x + "_") for x in contigs] 
        
        parallel_args = zip(new_pres, new_args)

        # star map will unpack the tuple and apply the args
        res = pool.starmap(func, parallel_args)

        # concat with pandas to drop header from each file
        df = pd.DataFrame()
        for fn in res:
            data = pd.read_table(fn,
                compression='gzip')
            df = pd.concat([df, data],axis=0,ignore_index=True)  
        
        output = outpre + "pileup_table.tsv.gz"
        df.to_csv(output, sep= "\t", index=False, header=True, compression='gzip')

        shutil.rmtree(tmp_dir, ignore_errors=False, onerror=None)
    
    # parse into bedgraphs (pos and neg)
    df = pd.read_table(output,
              compression='gzip')
    
    
    # parse into bedgraphs (pos and neg)
    
    df_pos = df[df.strand == "+"]
    df_neg = df[df.strand == "-"]

    format_bedgraphs(df_pos, depth, outpre + "pos_") 
    format_bedgraphs(df_neg, depth, outpre + "neg_")
    
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
    
    parser.add_argument('-o',
                        '--outpre',
                        help="""prefix for output files""",
                        required = False,
                        default = "")
    
    parser.add_argument('-t',
                        '--threads',
                        help="""Threads to use when running mpileup. If
                        threads is > 1 then the mpileup command will be
                        split up by chromosome to run using multiple
                        threads
                        default = 1""",
                        required = False,
                        default = 1,
                        type = int)
    
    parser.add_argument('-v',
                        '--verbose',
                        help="""print run information""",
                        required = False,
                        default = False)

    args=parser.parse_args()
    
    bam_name = args.bam

    # ok to have repeated args in samtools command
    pileup_args =  " -d 10000000 -L 10000000 -B " + args.pileup    
    fasta_name = args.fasta
    depth = args.depth
    outpre = args.outpre
    threads = args.threads
    verbose = args.verbose

    generate_mismatch_profile(bam_name, 
            fasta_name, 
            pileup_args,
            depth,
            outpre,
            threads,
            verbose) 
                
if __name__ == '__main__': main()

