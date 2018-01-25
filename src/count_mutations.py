#!/usr/bin/env python3
import argparse
import os 
import sys
import subprocess
import gzip
import pandas as pd
import uuid 
import shutil
import itertools
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

def generate_pileup(bam, fasta, min_depth, deletion_length, 
                    additional_args, outpre, region = "", verbose = False):
    """ returns fileobject to pileup output 
    
    args:
        bam = path to bam
        fasta = path to fasta
        additional_args = list of additional arguments for samtools pileup
    """

    output = open(outpre + "pileup_table.tsv.gz", "w")

    pileup_cmd =  "samtools view -h " + bam + " " + region + \
                  " | filterBam -d " + str(deletion_length) + \
                  " | samtools " + \
                  "mpileup " + \
                  "-f " + fasta + " " + \
                  additional_args + \
                  " - " + \
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

class Interval:
    """ bed interval object """
    def __init__(self, chrom, start, end, count):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.count = count
    def __str__(self):
        return "{}\t{}\t{}\t{}".format(self.chrom,
                self.start,
                self.end,
                self.count)

def convert_pileup(bg, outbg):
    """ convert pileup format to bedgraph 
    merging adjacent intervals with the same
    count values """ 
    
    ivl_cache = line_to_interval(bg.readline())
    current_chrom = ivl_cache.chrom

    for line in bg:
        
        ivl = line_to_interval(line)
        
        if ivl.chrom != current_chrom:
            
            outbg.write(str(ivl_cache) + "\n")
            ivl_cache = ivl
            current_chrom = ivl.chrom

        elif ivl_cache.end == ivl.start and ivl_cache.count == ivl.count:
            # same coverage, modify end
            ivl_cache.end = ivl.end
        
        else:
            outbg.write(str(ivl_cache) + "\n")
            ivl_cache = ivl

def line_to_interval(line):
    line = line.rstrip()
    fields = line.split("\t")
    return Interval(fields[0], 
                    int(fields[1]), 
                    int(fields[2]),
                    float(fields[3]))

def list_to_interval(lst):
    return Interval(lst[0], lst[1], lst[2], lst[3])

def merge_bedgraphs(prefix, strand, 
                    insuffix = ".bedgraph.tmp.gz",
                    outsuffix = ".bedgraph.gz"):
    def_fnames = [
      "mismatches",
      "insertions",
      "deletions",
      "depth"
    ]
    
    bgnames = [prefix + strand + x + insuffix for x in def_fnames]
    outnames = [prefix + strand + x + outsuffix for x in def_fnames]
    
    for idx, fn in enumerate(bgnames):
      outname = outnames[idx]
      with gzip.open(fn, 'rt') as f, gzip.open(outname, 'wt') as fout:
        convert_pileup(f, fout)
      os.unlink(fn)
    
def format_bedgraphs(df, depth, prefix):

    """ take pandas dataframe and generate bedgraphs """
    df = df.assign(mismatch_ratio = lambda df: df.mmcount / df.depth)
    df = df.assign(insertion_ratio = lambda df: df.inscount / df.depth)
    df = df.assign(deletion_ratio = lambda df: df.delcount / df.depth)

    df = df.assign(start = lambda df:df.pos - 1)
    df = df.rename(columns = {'pos':'end'})

    df = df[df.depth >= depth]

    df_mm = df[['chr', 'start', 'end', 'mismatch_ratio']]
    df_ins = df[['chr', 'start', 'end', 'insertion_ratio']]
    df_del = df[['chr', 'start', 'end', 'deletion_ratio']]
    df_depth = df[['chr', 'start', 'end', 'depth']]

    df_mm = df_mm.sort_values(['chr', 'start'], ascending=[True, True])
    df_ins = df_ins.sort_values(['chr', 'start'], ascending=[True, True])
    df_del = df_del.sort_values(['chr', 'start'], ascending=[True, True])
    df_depth = df_depth.sort_values(['chr', 'start'], ascending = [True, True])

    df_mm.to_csv(prefix + "mismatches.bedgraph.tmp.gz", sep = "\t", index=False, header=False, compression='gzip')
    df_ins.to_csv(prefix + "insertions.bedgraph.tmp.gz", sep = "\t", index=False, header=False, compression='gzip')
    df_del.to_csv(prefix + "deletions.bedgraph.tmp.gz", sep = "\t", index=False, header=False, compression='gzip')
    df_depth.to_csv(prefix + "depth.bedgraph.tmp.gz", sep = "\t", index=False, header=False, compression='gzip')

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
        threads, deletion_length, debug = False):
    
   if debug:
       print("started processing {}".format(str(datetime.now())),
           file = sys.stderr)
   
   outdir = os.path.dirname(outpre)
   if not os.path.exists(outdir) and outdir:
       os.makedirs(outdir)

   ## generate per nucleotide mismatch and indel counts
   if threads == 1:
       output = generate_pileup(input_bam, fasta, depth, 
                                deletion_length, additional_args,
                                outpre, region = "", verbose = debug)
   
   else:
       """ if multiple threads then run mpileup on each chromosome using
       the region arg passed to samtools view. 
       Write output to temp folder, then combine results. 
       """
       
       # generate list of new regional arguments to pass in parallel
       contigs = retrieve_header(input_bam)
       
       if "-r " in additional_args:
           print("-r option is not allowed when running with multiple threads", 
                 file = sys.stderr)
       
       

       pool = Pool(threads)
       
       # build function obj
       func = partial(generate_pileup,
           input_bam,
           fasta,
           depth,
           deletion_length,
           additional_args, 
           verbose = debug)
       
       # make tmp directory
       tmp_dir = "tmpfiles-" + str(uuid.uuid4())
       if not os.path.exists(tmp_dir): 
         os.makedirs(tmp_dir)
         print("temporary files placed in directory:\n{}".format(tmp_dir),
               file = sys.stdout)
         
       # generate prefixes that with tmp_dir and contig id
       new_pres = [os.path.join(tmp_dir, x + "_") for x in contigs] 
       
       parallel_args = zip(new_pres, contigs)
       
       ## star map will unpack the tuple and apply the args
       ## each produced filename will be returned in list
       res = []
       for results in pool.starmap(func, parallel_args):
           res.append(results)
       pool.close()
       pool.join()
       
       output = outpre + "pileup_table.tsv.gz"

       # combine gzipped per chromosome files
       output_tbl = gzip.open(output, 'wt')
       for idx, fn in enumerate(res):
           with gzip.open(fn, 'rt') as data:

               if idx == 0:
                   # keep header from first file
                   for line in data:
                       output_tbl.write(line)
               else:
                   # skip header for remaining
                   data.readline()
                   for line in data:
                       output_tbl.write(line)
       
       output_tbl.close()

       shutil.rmtree(tmp_dir, ignore_errors=False, onerror=None)
   
   # parse into bedgraphs (pos and neg)
   df = pd.read_table(output, compression = 'gzip')
   
   # parse into strand-specific bedgraphs (pos and neg)
   df_pos = df[df.strand == "+"]
   df_neg = df[df.strand == "-"]

   format_bedgraphs(df_pos, depth, outpre + "pos_") 
   format_bedgraphs(df_neg, depth, outpre + "neg_")
   
   merge_bedgraphs(outpre, "pos_")
   merge_bedgraphs(outpre, "neg_")
   
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
                        -L 10000000 (use up to 1e7 reads per base for indel calling)
                        do not count orphan reads (paired end)
                        do not double count if paired end reads overlap
                        """, 
                        required = False,
                        default = "")
    
    parser.add_argument('-d',
                        '--depth',
                        help = """minimum read coverage required for
                        reporting mismatch or indel frequencies. 
                        Default = 5""", 
                        required = False, 
                        default = 5, type = float)
                        
    parser.add_argument('-l',
                        '--deletion_length',
                        help = """do not count deletions greater than or 
                        equal to deletion_length, default = 4""",
                        default = 4, type = int)
                        
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

    args = parser.parse_args()
    
    bam_name = args.bam

    # ok to have repeated args in samtools command
    pileup_args =  " -d 10000000 -L 10000000 -B " + args.pileup    
    fasta_name = args.fasta
    depth = args.depth
    outpre = args.outpre
    deletion_length = args.deletion_length
    threads = args.threads
    verbose = args.verbose
    
    if not os.path.isfile(bam_name) or not os.path.isfile(fasta_name):
        sys.exit("input bam {} or fasta {} not found".format(bam_name,
            fasta_name))


    generate_mismatch_profile(bam_name, 
            fasta_name, 
            pileup_args,
            depth,
            outpre,
            threads,
            deletion_length,
            verbose) 
                
if __name__ == '__main__': main()

