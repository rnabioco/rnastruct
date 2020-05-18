#!/usr/bin/env python3
import argparse
import textwrap
import os 
import sys
import subprocess
import gzip
import pandas as pd
import numpy as np
import uuid 
import shutil
import itertools
import heapq
import io
import resource
import math
import multiprocessing as mp
import pysam
import atexit
import dask.dataframe as dd
import binascii
import pdb

from contextlib import ExitStack
from datetime import datetime
from collections import Counter
from functools import partial
from shutil import which
from operator import itemgetter

from utils import * 
import pileup_to_counts

bp_dict = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C"}



tbl_cols = [
 "chr",
 "pos",
 "strand", 
 "ref_base",
 "depth",
 "refcount",
 "acount",
 "ccount",
 "gcount",
 "tcount",
 "mmcount",
 "indelcount",
 "mismatch_ratio",
 "indel_ratio",
 "mutation_ratio",
 "stderr"]

def generate_pileup(bam, fasta, min_depth, deletion_length, 
                    additional_args, samflag, libtype, outpre, region = "", verbose = False):
    """ returns fileobject to pileup output 
    args:
        bam = path to bam
        fasta = path to fasta
        additional_args = list of additional arguments for mpileup
    """
    
    if libtype == "antisense": 
      rev_flag = " -r "
    else:
      rev_flag = ""
   
    if region is None:
        region = ""

    if not isinstance(region, list):
        region = [region]
        
    region = " ".join(region)
    
    output = outpre + "pileup_table.tsv.gz"
    
    output_tmpbam = outpre + "filteredbam.bam"

    bamfilter_cmd = "samtools view -h " + samflag + " " + bam + " " + region + \
                  " | filterBam -d " + str(deletion_length) + \
                  " | samtools view -b " + \
                  " > " + output_tmpbam + " ; " + \
                  " samtools index " + output_tmpbam

    filter_run = subprocess.run(bamfilter_cmd, 
            shell=True, 
            stderr = sys.stderr, 
            stdout = sys.stdout)

    if verbose:
        print("bamfilter command is:\n" + filter_run.args, 
                file = sys.stderr)
    mpileup_args = conv_args(additional_args)
#    pileup_args = []
#    for k,v in pileup_args.items():
#      pileup_args.append("{}={}".format(k,v))
#
#    pileup_args = " ".join(pileup_args)
#    pileup_cmd =  "bcftools " + \
#                  "mpileup " + \
#                  " -a AD " + \
#                  "-f " + fasta + " " + \
#                  mpileup_args + \
#                  " " + output_tmpbam + " " \
#                  " | " + \
#                  "pileup_to_counts.py -v -  -d " + \
#                  str(min_depth) + rev_flag + " -o " + output + \
#                  " -b " + output_tmpbam
#                  " --pileup-args " + pileup_args + " " 
    mpileup_args = mpileup_args.split()
    pileup_cmd =  ["bcftools", "mpileup", "-a", "AD", "-f", fasta] + \
                  mpileup_args + [output_tmpbam] 
    pileup_run = subprocess.Popen(pileup_cmd, 
            shell=False, 
            stderr = sys.stderr, 
            stdout = subprocess.PIPE)
    
    # cyvcf2 can use a file descriptor to open vcf 
    pileup_to_counts.vcf_to_counts(pileup_run.stdout.fileno(),
                                   output_tmpbam,
                                   output,
                                   additional_args,
                                   min_depth = min_depth,
                                   return_comp = libtype == "antisense",
                                   debug = False)
    
    if verbose:

        print("pileup command is:\n" + " ".join(pileup_run.args), 
                file =sys.stderr)
        print("formatted pileup output is here:\n" + output, 
                file = sys.stderr)
        print("completed pileup {}".format(str(datetime.now())), 
                file = sys.stderr)
    
    return output

class Interval:
    """ bed interval object """
    def __init__(self, chrom, start, end, values):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.vals = values
    def __str__(self):
        out_vals = [str(x) for x in self.vals]
        return "{}\t{}\t{}\t{}".format(self.chrom,
                self.start,
                self.end,
                "\t".join(out_vals))

def line_to_interval(line):
    line = line.rstrip()
    fields = line.split("\t")
    ivl = Interval(fields[0], 
                    int(fields[1]), 
                    int(fields[2]),
                    [float(x) for x in fields[3:]])
    return ivl

class Pileup:
    """ pileup object """
    def __init__(self, chrom, pos, strand, ref_base, aux):
        self.chrom = chrom
        self.pos = pos
        self.strand = strand
        self.ref_base = ref_base
        self.aux_fields = aux
        
    def __str__(self):
        aux_out = [str(x) for x in self.aux_fields]
        return "{}\t{}\t{}\t{}\t{}\t{}".format(self.chrom,
                self.pos,
                self.strand,
                aux_out[0],
                self.ref_base,
                "\t".join(aux_out[1:]))
                
    def ivl_index(self):
        return (self.chrom, self.pos, self.strand)

def line_to_pileup(line):
    line = line.rstrip()
    fields = line.split("\t")
    aux = [fields[3]] + fields[5:]
    aux = [int(x) for x in aux]
    return Pileup(fields[0], 
                    int(fields[1]), 
                    fields[2],
                    fields[4],
                    aux)

def file_to_pileup(fn):
    """ generator for pileup object from file """
    for line in fn:
        yield line_to_pileup(line)
        
def convert_pileup(bg, outbg):
    """ convert pileup format to bedgraph 
    merging adjacent intervals with the same
    count values
    """
    
    ivl_cache = line_to_interval(bg.readline())
    current_chrom = ivl_cache.chrom
    
    for line in bg:
        ivl = line_to_interval(line)
        if ivl.chrom != current_chrom:
            outbg.write("{}\n".format(str(ivl_cache)))
            ivl_cache = ivl
            current_chrom = ivl.chrom
        elif ivl_cache.end == ivl.start and ivl_cache.vals == ivl.vals:
            # same coverage, modify end
            ivl_cache.end = ivl.end
        else:
            outbg.write("{}\n".format(str(ivl_cache)))
            ivl_cache = ivl
            
    ## clear out the last interval        
    outbg.write("{}\n".format(str(ivl_cache)))
    

def process_bedgraph(fn, outname):

  if gz_is_empty(fn):
    print("{} is empty, no-data".format(fn), file = sys.stderr)
    if os.path.isfile(fn): os.unlink(fn)
    return     
  
  with gzip.open(outname, 'wt', compresslevel = 6) as fout:
    memmap_mm = io.StringIO()
    df = pd.read_csv(fn, sep = "\t", header = None)
    df = df.sort_values([0, 1, 2])
    df.to_csv(memmap_mm, sep = "\t", index = False, header = False)
    memmap_mm.seek(0)
    convert_pileup(memmap_mm, fout)
    memmap_mm.close()
  os.unlink(fn)

def merge_bedgraphs(prefix, strands, output_dir, threads = 1,  
                    insuffix = ".bedgraph.tmp",
                    outsuffix = ".bedgraph.gz"):
    def_fnames = [
      "mismatches",
      "indels",
      "depth",
      "mutations"
    ]
    
    bgnames = []
    outnames = []
    for strand in strands:
       bgnames += [os.path.join(prefix, strand + x + insuffix) for x in def_fnames]
       outnames += [output_dir + strand + x + outsuffix for x in def_fnames]
    
    pool = mp.Pool(threads)

    p = pool.starmap(process_bedgraph, zip(bgnames, outnames))
    pool.close()
    pool.join()

def write_bedgraphs(lst_dfs, lst_fns):
    """ write stringIO bedgraphs to supplied filenames"""
    
    for bg in zip(lst_dfs, lst_fns):
      with open(bg[1], 'a') as fd:
        bg[0].seek(0)
        fd.write(bg[0].getvalue())
        bg[0].close()

def memmap_df(df):    
    memmap_mm = io.StringIO()
    df.to_csv(memmap_mm, sep = "\t", index = False, header = False)
    memmap_mm.seek(0)
    df = io.StringIO()
    convert_pileup(memmap_mm, df)
    memmap_mm.close()
    return(df)
    
def format_bedgraphs(df, depth, nucs_to_keep, prefix):
    
    """ take pandas dataframe and split into bedgraph dataframes
    return list of dataframes and list of output filenames
    pass to write bedgraphs.
    
    """
    df = df.assign(mismatch_ratio = lambda df: df.mmcount / df.depth)
    df = df.assign(indel_ratio = lambda df: df.indelcount / df.depth)
    df = df.assign(mutation_ratio = lambda df: (df.mmcount + df.indelcount) / df.depth)
    df = df.assign(stderr = lambda df:
            (np.sqrt(df.mutation_ratio) / np.sqrt(df.depth)))

    df = df.assign(start = lambda df:df.pos - 1)
    df = df.rename(columns = {'pos':'end'})

    df = df[df.depth >= depth]
    df_depth = df[['chr', 'start', 'end', 'depth']]
    
    df = df[df.ref_base.isin(nucs_to_keep)]
    
    df_mm = df[['chr', 'start', 'end', 'mismatch_ratio']]
    df_indel = df[['chr', 'start', 'end', 'indel_ratio']]
    df_mut = df[['chr', 'start', 'end', 'mutation_ratio', 'stderr']]

    df_mm = df_mm.sort_values(['chr', 'start'], ascending=[True, True])
    df_indel = df_indel.sort_values(['chr', 'start'], ascending=[True, True])
    df_depth = df_depth.sort_values(['chr', 'start'], ascending = [True, True])
    df_mut = df_mut.sort_values(['chr', 'start'], ascending = [True, True])

    ## run merge intervals on memory mapped data
    df_list = [df_mm, df_indel, df_depth, df_mut]
    df_fnames = [
      "mismatches",
      "indels",
      "depth",
      "mutations"
    ]
     
    mem_mapped_dfs = []
    mem_mapped_dfs_fnames = []
    for idx, d in enumerate(df_list):
        if d.empty:
          continue
        
        ## convert to in-memory file and merge intervals in place
        fobj = memmap_df(d)
        mem_mapped_dfs.append(fobj)
        mem_mapped_dfs_fnames.append(df_fnames[idx])
    
    out_fns = [prefix + x + ".bedgraph.tmp" for x in mem_mapped_dfs_fnames]
    
    return(mem_mapped_dfs, out_fns, [df], [prefix + "pileup_all.tmp"])

def parse_library_type(bam_path, strandedness, libtype, skip_single_ended
        = False):
    """ return a list of appropriate flags for filtering bamfile 
      returns a list with a list for each bam flag setting
      first element is a filtering flag
      second element is a string that indicates if the alignments are sense or antisense
      
      if skip_single_ended, then will not count single-end alignments found
      mixed with paired-end alignments 
    """
    flags = {}
    if strandedness == "fr-firststrand":
      if libtype == "paired":
        flags["sense"] = []
        flags["antisense"] = []
        
        ## parse out forward strand alignments (R2) (second in pair, forward mapped)
        sam_flag_1 = "-f 128 -F 16 "
        flags["sense"].append(sam_flag_1)

        ## parse out forward strand alignments (R1) (first in pair reverse mapped (16 + 64))
        sam_flag_2 = "-f 80 "
        flags["sense"].append(sam_flag_2)
        
        ## parse out reverse alignments (R2) (second in pair, reverse mapped (128 + 16))
        sam_flag_3 = "-f 144 "
        flags["antisense"].append(sam_flag_3)

        ## parse out reverse alignments (R1) (first in pair, forward mapped)
        sam_flag_4 = "-f 64 -F 16"
        flags["antisense"].append(sam_flag_4)
        
        #output = [[sam_flag_1, lib_type_1], 
        #          [sam_flag_2, lib_type_2],
        #          [sam_flag_3, lib_type_3],
        #          [sam_flag_4, lib_type_4]]
        
        if not skip_single_ended:
            #output.append([ "-f 16 -F 1 ", "sense"])
            #output.append(["-F17 ", "antisense"])
           
            flags["sense"].append("-f 16 -F 1")
            flags["antisense"].append("-F17 ")

      else:
        # single end
        #output = []
        
        flags["sense"].append("-f 16 -F 1")
        flags["antisense"].append("-F17 ")
        
        #output.append([ "-f 16 -F 1 ", "sense"])
        #output.append(["-F17 ", "antisense"])
        
        
    elif strandedness == "fr-secondstrand":
      flags["sense"] = []
      flags["antisense"] = []
      
      if libtype == "paired":
        ## parse out forward strand alignments (R2) (second in pair, reverse mapped)
        sam_flag_1 = "-f 144 "
        flags["sense"].append(sam_flag_1)

        ## parse out forward strand alignments (R1) (first in pair not reverse mapped (16 + 64))
        sam_flag_2 = "-f 64 -F 16 "
        flags["sense"].append(sam_flag_2)
        
        ## parse out reverse alignments (R2) (second in pair, forward mapped)
        sam_flag_3 = "-f 128 -F 16 "
        flags["antisense"].append(sam_flag_3)
        
        ## parse out reverse alignments (R1) (first in pair, reverse mapped)
        sam_flag_4 = "-f 80"
        flags["antisense"].append(sam_flag_4)
        
        #output = [[sam_flag_1, lib_type_1], 
        #          [sam_flag_2, lib_type_2],
        #          [sam_flag_3, lib_type_3],
        #          [sam_flag_4, lib_type_4]]
        
        if not skip_single_ended:
            flags["antisense"].append("-f 16 -F 1")
            flags["sense"].append("-F17 ")
            #output.append([ "-f 16 -F 1 ", "antisense"])
            #output.append(["-F17 ", "sense"])

      else:
        # single end
        #output = []
        #output.append([ "-f 16 -F 1 ", "antisense"])
        #output.append(["-F17 ", "sense"])
        flags["antisense"].append("-f 16 -F 1")
        flags["sense"].append("-F17 ")
        
    elif strandedness == "unstranded":
      flags["sense"] = [""]
      
    else:
      sys.exit("libtype specification failed")
    
    return flags

def setdiff(lst_a, lst_b):
    sety = set(lst_b)
    return [x for x in lst_a if x not in sety]

def generate_mismatch_profile(input_bam, fasta, additional_args, depth, outpre, 
        threads, deletion_length, bam_flag, libtype, chroms_to_exclude = None,
        region_to_pileup = None, debug = False):
    
   if debug:
       print("started processing {}".format(str(datetime.now())),
           file = sys.stderr)

   chroms = retrieve_contigs(input_bam, require_reads = True)
   if chroms_to_exclude is not None:
     chroms = setdiff(chroms, chroms_to_exclude)

   ## generate per nucleotide mismatch and indel counts
   if threads == 1:
       if chroms_to_exclude is not None:
           region_to_pileup = chroms
       output = generate_pileup(input_bam, fasta, depth, 
                                deletion_length, additional_args,
                                bam_flag, libtype,
                                outpre, region = region_to_pileup, verbose = debug)
   
   else:
       """ if multiple threads then run mpileup on each chromosome using
       the region arg passed to samtools view. 
       Write output to temp folder, then combine results. 
       """
       
       # generate list of new regional arguments to pass in parallel
       #contigs = retrieve_header(input_bam)
       
       if region_to_pileup is not None:
           print("-r option is not respected when running with multiple threads", 
                 file = sys.stderr)
       
       pool = mp.Pool(threads)
       
       # build function obj
       func = partial(generate_pileup,
           input_bam,
           fasta,
           depth,
           deletion_length,
           additional_args, 
           bam_flag,
           libtype,
           verbose = debug)
        
       # generate prefixes that begin with tmp_dir and contig id
       new_pres = [os.path.join(outpre, x + "_") for x in chroms] 
       
       parallel_args = zip(new_pres, chroms)
       
       ## star map will unpack the tuple and apply the args
       ## each produced filename will be returned in list
       res = []
       for results in pool.starmap(func, parallel_args):
           res.append(results)
       pool.close()
       pool.join()
       
       output = outpre + "pileup_table.tsv.gz"

       # combine gzipped per chromosome files
       output_tbl = gzip.open(output, 'wt', compresslevel = 6)
       for idx, fn in enumerate(res):
           with gzip.open(fn, 'rt') as data:

               if idx == 0:
                   # keep header from first file
                   shutil.copyfileobj(data, output_tbl)
               else:
                   # skip header for remaining
                   data.readline()
                   shutil.copyfileobj(data, output_tbl)
       
       output_tbl.close()
   
   return(output)

def split_and_apply(df, min_depth, nucs_to_keep, outprefix):
    """ parse into strand-specific bedgraphs (pos and neg)
    and return a list of bedgraph dataframes and their output filenames.
    pass to write_bedgraph to write to disk
    """
    
    df_pos = df[df.strand == "+"]
    df_neg = df[df.strand == "-"]
    
    bgs_pos, bgs_pos_fns, pos_pileup, pileup_fn_pos = format_bedgraphs(df_pos, min_depth, 
            nucs_to_keep, os.path.join(outprefix, "pos_")) 
    bgs_neg, bgs_neg_fns, neg_pileup, pileup_fn_neg = format_bedgraphs(df_neg, min_depth, 
            nucs_to_keep, os.path.join(outprefix, "neg_"))
    
    bgs_out = bgs_pos + bgs_neg
    bgs_fns = bgs_pos_fns + bgs_neg_fns
    pileup_dfs = pos_pileup + neg_pileup
    pileup_fn = pileup_fn_pos + pileup_fn_neg
    return([bgs_out, bgs_fns, pileup_dfs, pileup_fn])
    
def write_pileup(df, output_fn):
    
    for idx,d in enumerate(df):
      d.to_csv(output_fn[idx], 
               mode = 'a', 
               sep = "\t",
               header = False, 
               index = False)

def write_pileup_df(df, output_fn):
    
    df.to_csv(output_fn, 
              compression = "gzip",
              mode = 'a', 
              sep = "\t",
              header = False, 
              index = False)

def generate_bedgraphs(pileup_fn, depth, outpre, threads, nucleotides,
        chunk_size = 1000000):
   """ master function for generating bedgraphs in parallel """
   
   ## read in pileup table in chunks to keep memory low
   reader = pd.read_hdf(pileup_fn, chunksize = chunk_size, iterator = True)
   
   ## apply multithreading if cpus available
   pool = mp.Pool(threads)
   
   func = partial(split_and_apply,
            min_depth = depth,
            nucs_to_keep = nucleotides,
            outprefix = outpre)
  
   ## results is a list of pandas df's and output filenames
   
   ## might need to check and remove preexisting res[3] filenames

   for res in pool.imap(func, reader):
       write_bedgraphs(res[0], res[1])
       write_pileup(res[2], res[3])

def mismatch_stats(df, min_depth, nucs_to_keep):
    
    df = df.assign(mismatch_ratio = lambda df: df.mmcount / df.depth)
    df = df.assign(indel_ratio = lambda df: df.indelcount / df.depth)
    df = df.assign(mutation_ratio = lambda df: (df.mmcount + df.indelcount) / df.depth)
    df = df.assign(stderr = lambda df:
            (np.sqrt(df.mutation_ratio) / np.sqrt(df.depth)))
    
    df = df[df.depth >= min_depth]
    df = df[df.ref_base.isin(nucs_to_keep)]
    
    return df

def calc_mismatch_stats(pileup_fn, depth, outpre, threads, nucleotides,
        chunk_size = 1000000):
   
   output_fn = os.path.splitext(pileup_fn)[0] + ".tmp"
   
   ## read in pileup table in chunks to keep memory low
   try:
       reader = pd.read_csv(pileup_fn, 
           sep = "\t",
           compression = 'gzip',
           chunksize = chunk_size, 
           iterator = True)
   except pd.errors.EmptyDataError:
       fout = gzip.open(output_fn, 'wt')
       fout.close()
       return output_fn

   ## apply multithreading if cpus available
   pool = mp.Pool(threads)
   
   func = partial(mismatch_stats,
            min_depth = depth,
            nucs_to_keep = nucleotides)
   
   ## res is a list of pandas df's and output filenames
   for res in pool.imap(func, reader):
       write_pileup_df(res, output_fn)
   
   return output_fn

def compute_chunks(n_lines, n_default_lines = 100000):
    """
    check ulimit and split lines into 
    proper number of chunks without exceed ulimit
    
    make each chunk n_default_lines unless 90% of ulimit would be
    exceeded
    """

    softlimit, hardlimit = resource.getrlimit(resource.RLIMIT_NOFILE)

    n_chunk_lines = int(n_lines / int(softlimit * 0.9)) 
    if n_chunk_lines < n_default_lines: 
        n_chunk_lines = n_default_lines
    
    n_chunks = math.ceil(n_lines / n_chunk_lines)

    return n_chunk_lines, n_chunks


def merge_pileup_tables(pileup_fns, output_pileup_fn, tmp_dir, 
        threads = 1, out_colnames = tbl_cols, verbose = True):
    """ sort and merge the sense and antisense pileup tables
    samtools returns pileup format sorted by largest chromosome in 
    descending order
    """
    
    if verbose:
        print("started concatenating anti-sense and sense pileup tables",
              file = sys.stderr)
        
    ## concat the pileup tables
    output_tmp_fn = os.path.join(tmp_dir, "tmp_table.tsv")
    with gzip.open(output_tmp_fn, 'wt', compresslevel = 6) as out_fh:
      out_fh.write("\t".join(out_colnames) + "\n")
      for idx, pfn in enumerate(pileup_fns):
        fh = gzip.open(pfn, 'rt')
        shutil.copyfileobj(fh, out_fh)
        fh.close()
    
    ## sort 
    out_tmp_sorted_ptable = unix_sort(output_tmp_fn, 
            output_pileup_fn, 
            threads,
            memory = "8G", 
            verbose = verbose)

def format_tbls(fn1, fn2, outfn):
     
    df_pos = pd.read_csv(fn1, sep = "\t", header = None, names = tbl_cols)
    df_neg = pd.read_csv(fn2, sep = "\t", header = None, names = tbl_cols)

    df = pd.concat([df_pos, df_neg])
    df = df.drop('start', 1)
    df = df.sort_values(['chr', 'pos', 'strand'], ascending = [True, True, True])

    df.to_csv(outfn, index = False, sep = "\t", compression = "gzip")

def main():
    
    parser = argparse.ArgumentParser(description="""
    Parse bam file and enumerate mismatches, insertions, and deletions per
    nucleotide. Generates bedgraphs for mismatches, insertions and
    deletions, as well as a summary table with aggregated information""",
    formatter_class = argparse.RawTextHelpFormatter )

    parser.add_argument('-b',
                        '--bam',
                        help ="""indexed bam file input,
                        by default the bam will be split into
                        forward and reverse bams,
                        pass the forward and reverse bams
                        to bypass splitting the bams
                        (i.e. forward.bam,reverse.bam)
                        \n""",
                        required = True)
    parser.add_argument('-f',
                        '--fasta',
                        help = textwrap.dedent("""\
                        path to fasta file indexed with samtools faidx, 
                        passed to samtools mpileup
                        \n"""),
                        required = True)
                        
    parser.add_argument('-L',
                        '--library',
                        help = """library type: either 'paired' or 'single' (default: %(default)s)
                        \n""",
                        required = False,
                        default = 'paired')
                        
    parser.add_argument('-s',
                        '--strandedness',
                        help = textwrap.dedent("""\
                        strandedness of library:
                        
                          'fr-firststrand' = first strand sequenced as R1 (default)
                             i.e. R1 alignments are the reverse complement of RNA, 
                             for paired-end R2 alignments are the same sequence as the RNA
                             
                          'fr-secondstrand' = second strand sequenced as R1 
                             i.e. R1 alignments are the same sequence as the RNA
                             for paired-end R2 alignments are the reverse complement of the RNA
                             
                          'unstranded' = report strandedness without respect to R1 or R2
                          (default: %(default)s)
                          \n"""),
                        required = False,
                        default = 'fr-firststrand')
                        
    parser.add_argument('-d',
                        '--depth',
                        help = textwrap.dedent("""\
                        minimum read coverage required for
                        reporting mismatch or indel frequencies (default: %(default)s)
                        \n"""), 
                        required = False, 
                        default = 5, 
                        type = int)
                        
    parser.add_argument('-l',
                        '--deletion_length',
                        help = textwrap.dedent("""\
                        do not count deletions greater than or equal 
                        to deletion_length (default: %(default)s)
                        \n"""),
                        required = False,
                        default = 4, type = int)
                        
    parser.add_argument('-o',
                        '--outpre',
                        help="""prefix for output files
                        \n""",
                        required = False,
                        default = "")
    
    parser.add_argument('-t',
                        '--threads',
                        help=textwrap.dedent("""\
                        Threads to use when running mpileup. If
                        threads is > 1 then the mpileup command will be
                        split up by chromosome to run using multiple
                        threads (default: %(default)s)
                        \n"""),
                        required = False,
                        default = 1,
                        type = int)
                        
    parser.add_argument('-n',
                        '--nucleotides',
                        help=textwrap.dedent("""\
                        Nucleotides to use for computing
                        mismatch and indel ratios. Provide
                        as a string. i.e. to report 
                        for all nucleotides. "ATCG"
                        (default: %(default)s)
                        \n"""),
                        required = False,
                        default = "AC",
                        type = str)
                        
    parser.add_argument('-ss',
                        '--skip-single-reads',
                        help=textwrap.dedent("""\
                        If set, skip single end reads when encountered
                        in a paired-end library.
                        (default: %(default)s , single end reads are
                        counted by default)
                        \n"""),
                        action = 'store_true')
                        
    parser.add_argument('-c',
                        '--chroms-to-skip',
                        help=textwrap.dedent("""\
                        space separated list of chroms to ignore
                        (default: %(default)s)
                        \n"""),
                        required = False,
                        nargs = "+",
                        default = None)
    parser.add_argument('-r',
                        '--region',
                        help=textwrap.dedent("""\
                        region to query (as samtools region string)
                        (default: %(default)s)
                        \n"""),
                        required = False,
                        default = None)
    parser.add_argument('-v',
                        '--verbose',
                        help="""print run information (default: %(default)s)\n""",
                        action = 'store_true')
    
    parser.add_argument('-k',
                        '--keep-temp-files',
                        help="""don't delete temp files (default: %(default)s)\n""",
                        action = 'store_false')
    
    parser.add_argument("--pileup-args",
                      nargs='+',
                      action=kvdictAppendAction,
                      metavar="KEY=VALUE",
                      help=textwrap.dedent("""\
                      Add key/value pileup arguments to overwrite
                      defaults, also can use --pileup-arg-fn to specify args\n"""
                      ))
    parser.add_argument("--pileup-arg-fn",
                      help=textwrap.dedent("""\
                      specify custom pileup.yaml file to overwrite default
                      arguments\n"""
                      ),
                      required = False)
    
    args = parser.parse_args()
    
    bam_name = args.bam
    split_bams = True
    if len(bam_name.split(",")) == 2:
       split_bams = False

    ### check system exectuables
    
    if not is_tool("samtools"):
        sys.exit("samtools is not in path")
    if not is_tool("bcftools"):
        sys.exit("bcftools is not in path")
    if not is_tool("pileup_to_counts.py"):
        sys.exit("pileup_to_counts.py is not in path, please add rnastruct/src to your PATH")
    if not is_tool("filterBam"):
        sys.exit("filterBam is not in path, please add rnastruct/src to your PATH")

    
    pileup_args = get_pileup_args(args.pileup_arg_fn, args.pileup_args)
    fasta_name = args.fasta
    depth = args.depth
    outpre = args.outpre
    deletion_length = args.deletion_length
    threads = args.threads
    verbose = args.verbose
    library = args.library
    nucleotides = args.nucleotides
    strandedness = args.strandedness
    skip_singles = args.skip_single_reads
    chroms_to_skip = args.chroms_to_skip
    region = args.region
   
    #### check options
    if not os.path.isfile(fasta_name):
        sys.exit("input fasta {} not found".format(fasta_name))

    if not os.path.isfile(bam_name):
        found = False
        if not split_bams:
            found = all([os.path.isfile(x) for x in bam_name.split(",")])
        if not found:
            sys.exit("input bam {} not found".format(bam_name))
    
    library_opts = ['single', 'paired']
    stranded_opts = ['fr-firststrand', 'fr-secondstrand', 'unstranded']
    
    if library not in library_opts:
        mess = "unknown option for --library: {} \ntry one of: {}"
        sys.exit(mess.format(library, ",".join(library_opts)))
    
    if strandedness not in stranded_opts:
        mess = "unknown option for --strandedness: {} \ntry one of: {}"
        sys.exit(mess.format(strandedness, ",".join(stranded_opts)))
     
    nucleotides = [x.upper() for x in nucleotides]

    ## keep the thread count reasonable
    threads = min(mp.cpu_count(), threads)
    
    #### make directories
    outdir = os.path.dirname(outpre)
    if not os.path.exists(outdir) and outdir:
       os.makedirs(outdir)
       
    ## make tmp directory
    tmp_dir = os.path.join(outdir, "tmpfiles")
    if not os.path.exists(tmp_dir): 
      os.makedirs(tmp_dir)
      if verbose:
          print("temporary files placed in directory:\n{}".format(tmp_dir),
                file = sys.stderr)
    else:
      sys.exit("temporary files directory already exists:\n{}".format(tmp_dir))
    
    ## set cleanup functions
    atexit.register(cleanup, tmp_dir, delete = args.keep_temp_files)

    #### parse library type  
    if split_bams:
        bam_flags = parse_library_type(bam_name, strandedness, library,
            skip_singles)
    else:
        bam_flags = {}
        bam_flags["sense"] = bam_name.split(",")[0]
        bam_flags["antisense"] = bam_name.split(",")[1]

    ## pass bams with proper flags to filter to pos and neg alignments
    ## based on library type
    ## merge tables and filter for depth for sense and antisense
    ## independently
    pileup_tbls = []
    
    for align_type,flags in bam_flags.items():
      tmp_pileup_tbls = []
      
      # split out antisense/sense bam as tmp file
      if split_bams:
          tmp_bam = os.path.join(tmp_dir, align_type + ".bam")
          split_bam(bam_name, tmp_bam, flags, threads = threads, memory = "1G")
      else:
          tmp_bam = flags
      
      pileup_fn = os.path.join(tmp_dir, align_type)
      if not os.path.exists(pileup_fn): 
        os.makedirs(pileup_fn)
          
      pileup_tbl = generate_mismatch_profile(tmp_bam, 
          fasta_name, 
          pileup_args,
          depth,
          pileup_fn,
          threads,
          deletion_length,
          " ",
          align_type,
          chroms_to_skip,
          region,
          verbose) 
          
      pileup_tbls.append(pileup_tbl)

    
    new_pileups = []
    for pileup in pileup_tbls:
      new_pileup = calc_mismatch_stats(pileup, 
                          depth,
                          tmp_dir,
                          threads,
                          nucleotides)
      new_pileups.append(new_pileup)

    #pdb.set_trace()

    print("merging and sorting pileup tables", file = sys.stderr)
    output_pileup_fn = outpre + "pileup_table.tsv.gz"
    merge_pileup_tables(new_pileups, output_pileup_fn, tmp_dir, 
            threads, verbose = verbose)

    #### bgzip and index
    out_bgzip = bgzip(output_pileup_fn)        
    pysam.tabix_index(out_bgzip, 
                      seq_col = 0, 
                      start_col = 1,
                      end_col = 1, 
                      zerobased = False,
                      force = True, 
                      line_skip = 1)

    print("Done", file = sys.stderr)
                
if __name__ == '__main__': main()

