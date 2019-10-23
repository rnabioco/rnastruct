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

from contextlib import ExitStack
from datetime import datetime
from collections import Counter
from functools import partial
from shutil import which
from operator import itemgetter

from utils import * 


bp_dict = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C"}


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
    

    pileup_cmd =  "samtools view -h " + samflag + " " + bam + " " + region + \
                  " | filterBam -d " + str(deletion_length) + \
                  " | bcftools " + \
                  "mpileup " + \
                  "-f " + fasta + " " + \
                  additional_args + \
                  " - " + \
                  " | " + \
                  "pileup_to_counts.py -v -  -d " + \
                  str(min_depth) + rev_flag + " -o " + output
                 
    pileup_run = subprocess.run(pileup_cmd, 
            shell=True, 
            stderr = sys.stderr, 
            stdout = sys.stdout)
    
    if verbose:

        print("pileup command is:\n" + pileup_run.args, 
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
    
    if strandedness == "fr-firststrand":
      if libtype == "paired":
        ## parse out forward strand alignments (R2) (second in pair, forward mapped)
        sam_flag_1 = "-f 128 -F 16 "
        lib_type_1 = "sense"
        ## parse out forward strand alignments (R1) (first in pair reverse mapped (16 + 64))
        sam_flag_2 = "-f 80 "
        lib_type_2 = "sense"
        
        ## parse out reverse alignments (R2) (second in pair, reverse mapped (128 + 16))
        sam_flag_3 = "-f 144 "
        lib_type_3 = "antisense"
        ## parse out reverse alignments (R1) (first in pair, forward mapped)
        sam_flag_4 = "-f 64 -F 16"
        lib_type_4 = "antisense"
        
        output = [[sam_flag_1, lib_type_1], 
                  [sam_flag_2, lib_type_2],
                  [sam_flag_3, lib_type_3],
                  [sam_flag_4, lib_type_4]]
        
        if not skip_single_ended:
            output.append([ "-f 16 -F 1 ", "sense"])
            output.append(["-F17 ", "antisense"])
        return output

      else:
        # single end
        output = []
        output.append([ "-f 16 -F 1 ", "sense"])
        output.append(["-F17 ", "antisense"])
        return output
        
        
    elif strandedness == "fr-secondstrand":
      if libtype == "paired":
        ## parse out forward strand alignments (R2) (second in pair, reverse mapped)
        sam_flag_1 = "-f 144 "
        lib_type_1 = "sense"
        ## parse out forward strand alignments (R1) (first in pair not reverse mapped (16 + 64))
        sam_flag_2 = "-f 64 -F 16 "
        lib_type_2 = "sense"
        
        ## parse out reverse alignments (R2) (second in pair, forward mapped)
        sam_flag_3 = "-f 128 -F 16 "
        lib_type_3 = "antisense"
        ## parse out reverse alignments (R1) (first in pair, reverse mapped)
        sam_flag_4 = "-f 80"
        lib_type_4 = "antisense"
        
        output = [[sam_flag_1, lib_type_1], 
                  [sam_flag_2, lib_type_2],
                  [sam_flag_3, lib_type_3],
                  [sam_flag_4, lib_type_4]]
        
        if not skip_single_ended:
            output.append([ "-f 16 -F 1 ", "antisense"])
            output.append(["-F17 ", "sense"])
        return output

      else:
        # single end
        output = []
        output.append([ "-f 16 -F 1 ", "antisense"])
        output.append(["-F17 ", "sense"])
        return output
        
    elif strandedness == "unstranded":
      
      return [["", "sense"]]
      
    else:
      sys.exit("libtype specification failed")

def setdiff(lst_a, lst_b):
    sety = set(lst_b)
    return [x for x in lst_a if x not in sety]

def generate_mismatch_profile(input_bam, fasta, additional_args, depth, outpre, 
        threads, deletion_length, bam_flag, libtype, chroms_to_exclude = None,
        debug = False):
    
   if debug:
       print("started processing {}".format(str(datetime.now())),
           file = sys.stderr)

   chroms = retrieve_contigs(input_bam, require_reads = True)
   if chroms_to_exclude is not None:
     chroms = setdiff(chroms, chroms_to_exclude)

   ## generate per nucleotide mismatch and indel counts
   if threads == 1:
       if "-r " in additional_args:
           ## pass region to samtools view
           args = additional_args.split()
           region_idx = args.index("-r")
           region_to_pileup = args[region_idx + 1]
           
           ## don't pass region through additional args
           ## samtools mpileup wont work with regional arg (input is sam)
           args.remove("-r")
           args.remove(region_to_pileup)
           additional_args = " ".join(args)
       else:
           if chroms_to_exclude is not None:
             region_to_pileup = chroms
           else:
             region_to_pileup = None
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
       
       if "-r " in additional_args:
           print("-r option is not allowed when running with multiple threads", 
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
        threads = 1, verbose = True):
    """ sort and merge the sense and antisense pileup tables
    samtools returns pileup format sorted by largest chromosome in 
    descending order
    values for duplicated intervals then can be summed in one pass
    see https://stackoverflow.com/questions/23450145/sort-a-big-file-with-python-heapq-merge
    """
    
    if verbose:
        print("started concatenating anti-sense and sense pileup tables",
              file = sys.stderr)
              
    ## concat the pileup tables
    out_tmp_ptable = os.path.join(tmp_dir, "tmp_pileup_table.tsv.gz")
    
    with gzip.open(out_tmp_ptable, 'wt', compresslevel = 6) as out_tmp_ptable_fn:
      for idx, pfn in enumerate(pileup_fns):
        fn = gzip.open(pfn, 'rt')
        
        if idx == 0:
          header = fn.readline()
          shutil.copyfileobj(fn, out_tmp_ptable_fn)
        else:
          tmp = fn.readline()
          shutil.copyfileobj(fn, out_tmp_ptable_fn)
        
        fn.close()
    
    out_tmp_sorted_ptable = unix_sort(out_tmp_ptable, tmp_dir, threads,
            memory = "8G", verbose = verbose)

    os.unlink(out_tmp_ptable) 
    if verbose:
        print("summing up columns from sorted tables",
              file = sys.stderr)
    
    cols = header.rstrip().split("\t")
    df = dd.read_csv(out_tmp_sorted_ptable, names = cols, sep = "\t")
    df = df.groupby(['chr', 'pos', 'strand', 'ref_base']).sum().compute(num_workers = threads)
    df = df.reset_index()
    df.to_hdf(output_pileup_fn, 'df', format = 'table')

    os.unlink(out_tmp_sorted_ptable)
    for f in pileup_fns:
      os.unlink(f)

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
 "stderr",
 "start"]

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
                        help ="""indexed bam file input
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
                        
    parser.add_argument('-p',
                        '--pileup',
                        help = textwrap.dedent("""\
                        additional command line arguments to pass to samtools mpileup 
                        by default -f is set by the --fasta argument to this script
                        
                        The following arguments are set by default, but can be modified.
                        --ff UNMAP,SECONDARY,QCFAIL,DUP (filter alignments)
                        -B (disable BAQ calculation) 
                        -d 1000000 (use up to 1e6 reads per base)
                        -L 1000000 (use up to 1e6 reads per base for indel calling)
                        -A count orphan reads (paired end)
                        -x disable read-pair overlap detection
                        \n"""), 
                        required = False,
                        default = "")
    
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

    parser.add_argument('-v',
                        '--verbose',
                        help="""print run information (default: %(default)s)\n""",
                        action = 'store_true')

    parser.add_argument('-k',
                        '--keep-temp-files',
                        help="""don't delete temp files (default: %(default)s)\n""",
                        action = 'store_false')
    
    parser.add_argument('--default-pileup-args',
                        help = textwrap.dedent("""\
                        The following arguments are set by default
                        --ff UNMAP,SECONDARY,QCFAIL,DUP (filter alignments)
                        -B (disable BAQ calculation)
                        -d 1000000 (use up to 1e6 reads per base)
                        -I dont call indels
                        -A count orphan reads (paired end)
                        -x disable read-pair overlap detection
                        -a AD
                        pass a string to replace these args
                        (default: %(default)s)\n"""
                        ),
                        required = False,
                        default = " --ff UNMAP,SECONDARY,QCFAIL,DUP -a AD -A -x -d 100000 -L 100000 -B -O v ")
    
    args = parser.parse_args()
    
    bam_name = args.bam
    
    ### check system exectuables
    
    if not is_tool("samtools"):
        sys.exit("samtools is not in path")
    if not is_tool("pileup_to_counts.py"):
        sys.exit("pileup_to_counts.py is not in path, please add rnastruct/src to your PATH")
    if not is_tool("filterBam"):
        sys.exit("filterBam is not in path, please add rnastruct/src to your PATH")

    # ok to have repeated args in samtools command
    pileup_args = args.default_pileup_args + " " +args.pileup    
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

    #### check options
    if not os.path.isfile(bam_name) or not os.path.isfile(fasta_name):
        sys.exit("input bam {} or fasta {} not found".format(bam_name,
            fasta_name))
    
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
    tmp_dir = os.path.join(outdir, "tmpfiles-" + str(uuid.uuid4()))
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
    bam_flags = parse_library_type(bam_name, strandedness, library,
            skip_singles)
    
    ## parse bam into two new bams if paired end
    ## one bam with all alignments that report the correct strand of the fragment
    ## the other bam with all alignments that are rev-comp of the fragment
    ## invert the second bam reported strands and merge
    pileup_tbls = []
      
    for idx, bam_flag in enumerate(bam_flags):
      bam_flag_filter = bam_flag[0]
      align_type = bam_flag[1]
      pileup_fn = os.path.join(tmp_dir, bam_flag[1]+ "_" + str(idx))
      if not os.path.exists(pileup_fn): 
        os.makedirs(pileup_fn)
        
      pileup_tbl = generate_mismatch_profile(bam_name, 
        fasta_name, 
        pileup_args,
        depth,
        pileup_fn,
        threads,
        deletion_length,
        bam_flag_filter,
        align_type,
        chroms_to_skip,
        verbose) 
      
      pileup_tbls.append(pileup_tbl)
        
    output_pileup_fn = os.path.join(tmp_dir, "pileup_table.hd5")
    merge_pileup_tables(pileup_tbls, output_pileup_fn, tmp_dir, threads, verbose)
      
    
    print("parsing pileup into bedgraph format", file = sys.stderr)

    ## parse output into bedgraphs
    generate_bedgraphs(output_pileup_fn, depth, tmp_dir, threads, nucleotides)
    
    format_tbls(os.path.join(tmp_dir, "pos_pileup_all.tmp"), 
                os.path.join(tmp_dir, "neg_pileup_all.tmp"),
                outpre + "pileup_table.tsv.gz")

    print("merging redundant bedgraph entries", file = sys.stderr)
    
    ## merge redundant intervals
    merge_bedgraphs(tmp_dir, ["pos_", "neg_"], outpre, threads)
    
    ## bgzip and index
    out_tbl = outpre + "pileup_table.tsv.gz" 
    out_bgzip = bgzip(out_tbl)        
    pysam.tabix_index(out_bgzip, 
                      seq_col = 0, 
                      start_col = 1,
                      end_col = 1, 
                      zerobased = False,
                      force = True, 
                      line_skip = 1)

    print("Done", file = sys.stderr)
                
if __name__ == '__main__': main()

