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

from contextlib import ExitStack
from datetime import datetime
from collections import Counter
from functools import partial
from shutil import which
from operator import itemgetter

def cleanup(tmp_dir, delete = True):
    
    if delete:
      print("removing temp directory: {}".format(tmp_dir))
      shutil.rmtree(tmp_dir, ignore_errors=False, onerror=None)     
    
def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    return which(name) is not None

bp_dict = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C"}

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
                    additional_args, samflag, libtype, outpre, region = "", verbose = False):
    """ returns fileobject to pileup output 
    args:
        bam = path to bam
        fasta = path to fasta
        additional_args = list of additional arguments for samtools pileup
    """
    
    if libtype == "antisense": 
      rev_flag = " -r"
    else:
      rev_flag = ""
   
    if region is None:
        region = ""

    if not isinstance(region, list):
       # tmp_region = []
       # tmp_region.append(region)
        region = [region]
        
    # quote to protect strange chrom names (i.e. rRNA gi|555853|gb|U13369.1|HSU13369) 
    for idx,val in enumerate(region):
      if not val.startswith("'") and val != "":
        region[idx] = "'" + val + "'"
    
    region = " ".join(region)
    output = open(outpre + "pileup_table.tsv.gz", "w")
    

    pileup_cmd =  "samtools view -h " + samflag + " " + bam + " " + region + \
                  " | filterBam -d " + str(deletion_length) + \
                  " | samtools " + \
                  "mpileup " + \
                  "-f " + fasta + " " + \
                  additional_args + \
                  " - " + \
                  " | " + \
                  "mpileupToReadCounts -d " + \
                  str(min_depth) + rev_flag + \
                  " | gzip "
                 
    if verbose:

        print("pileup command is:\n" + pileup_cmd, file = sys.stderr)
        
        pileup_run = subprocess.run(pileup_cmd, 
            shell=True, stderr = sys.stderr, stdout = output)
        
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
    
def gz_is_empty(fname):
    ''' Test if gzip file fname is empty
        Return True if the uncompressed data in fname has zero length
        or if fname itself has zero length
        Raises OSError if fname has non-zero length and is not a gzip file
        https://stackoverflow.com/questions/37874936/how-to-check-empty-gzip-file-in-python
    '''
    if not os.path.isfile(fname) :
      return True
    
    if fname.endswith(".gz"):
        with gzip.open(fname, 'rb') as f:
            data = f.read(1)
        return len(data) == 0
    else: 
       with open(fname, 'r') as f:
            data = f.read(1)
       return len(data) == 0

def merge_bedgraphs(prefix, strand, output_dir, 
                    insuffix = ".bedgraph.tmp",
                    outsuffix = ".bedgraph.gz"):
    def_fnames = [
      "mismatches",
      "insertions",
      "deletions",
      "depth",
      "mutations"
    ]
    
    bgnames = [os.path.join(prefix, strand + x + insuffix) for x in def_fnames]
    outnames = [output_dir + strand + x + outsuffix for x in def_fnames]
    
    for idx, fn in enumerate(bgnames):
      outname = outnames[idx]
      if gz_is_empty(fn):
        print("{} is empty, no-data".format(fn), file = sys.stderr)
        if os.path.isfile(fn): os.unlink(fn)
        continue
      
      with gzip.open(outname, 'wt') as fout:
        memmap_mm = io.StringIO()
        df = pd.read_table(fn, header = None)
        df = df.sort_values([0, 1, 2])
        df.to_csv(memmap_mm, sep = "\t", index = False, header = False)
        memmap_mm.seek(0)
        convert_pileup(memmap_mm, fout)
        memmap_mm.close()
      os.unlink(fn)

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
    df = df.assign(insertion_ratio = lambda df: df.inscount / df.depth)
    df = df.assign(deletion_ratio = lambda df: df.delcount / df.depth)
    df = df.assign(mutation_ratio = lambda df: (df.mmcount + df.delcount) / df.depth)
    df = df.assign(stderr = lambda df:
            (np.sqrt(df.mutation_ratio) / np.sqrt(df.depth)))

    df = df.assign(start = lambda df:df.pos - 1)
    df = df.rename(columns = {'pos':'end'})

    df = df[df.depth >= depth]
    df_depth = df[['chr', 'start', 'end', 'depth']]
    
    df = df[df.ref_base.isin(nucs_to_keep)]
    
    df_mm = df[['chr', 'start', 'end', 'mismatch_ratio']]
    df_ins = df[['chr', 'start', 'end', 'insertion_ratio']]
    df_del = df[['chr', 'start', 'end', 'deletion_ratio']]
    df_mut = df[['chr', 'start', 'end', 'mutation_ratio', 'stderr']]

    df_mm = df_mm.sort_values(['chr', 'start'], ascending=[True, True])
    df_ins = df_ins.sort_values(['chr', 'start'], ascending=[True, True])
    df_del = df_del.sort_values(['chr', 'start'], ascending=[True, True])
    df_depth = df_depth.sort_values(['chr', 'start'], ascending = [True, True])
    df_mut = df_mut.sort_values(['chr', 'start'], ascending = [True, True])

    ## run merge intervals on memory mapped data
    df_list = [df_mm, df_ins, df_del, df_depth, df_mut]
    df_fnames = [
      "mismatches",
      "insertions",
      "deletions",
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
        ## parse out sense alignments (R2) (second in pair)
        sam_flag_1 = "-f 128 "
        lib_type_1 = "sense"
        ## parse out antisense alignments (R1) (first in pair)
        sam_flag_2 = "-f 64 "
        lib_type_2 = "antisense"
        
        output = [[sam_flag_1, lib_type_1], [sam_flag_2, lib_type_2]]
        
        if not skip_single_ended:
            output.append(["-F 1", "antisense"])
        return output

      else:
        # single end
        return [["", "antisense"]]
        
    elif strandedness == "fr-secondstrand":
      if libtype == "paired":
        ## parse out sense alignments (R1) (firsti n pair)
        sam_flag_1 = "-f 64 "
        lib_type_1 = "sense"
        ## parse out antisense alignments (R2) (second in pair)
        sam_flag_2 = "-f 128 "
        lib_type_2 = "antisense"
        
        output = [[sam_flag_1, lib_type_1], [sam_flag_2, lib_type_2]]
        
        if not skip_single_ended:
            output.append(["-F 1", "sense"])
        return output

      else:
        # single end
        return [["", "sense"]] 
        
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

   chroms = retrieve_header(input_bam)
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
      d.to_csv(output_fn[idx], mode = 'a', sep = "\t",
              index = False)

def generate_bedgraphs(pileup_fn, depth, outpre, threads, nucleotides, chunk_size = 100000):
   """ master function for generating bedgraphs in parallel """
   
   ## read in pileup table in chunks to keep memory low
   reader = pd.read_table(pileup_fn, compression = 'gzip', chunksize = chunk_size)
   
   ## apply multithreading if cpus available
   pool = mp.Pool(threads)
   
   func = partial(split_and_apply,
            min_depth = depth,
            nucs_to_keep = nucleotides,
            outprefix = outpre)
  
   ## results is a list of pandas df's and output filenames
   results = pool.imap(func, reader)
   pool.close()
   pool.join()
   del reader
   
   ## need to remove preexisting res[3] filenames

   for res in results:
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

def ungzip(in_fn, out_fn):

    with gzip.open(in_fn, 'rb') as f_in:
      with open(out_fn, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    
    os.unlink(in_fn)

def bgzip(in_fn):
    """
    convert file to bgzipped format
    """
    tmp_out_fn = in_fn.replace(".gz", "")
    out_fn = in_fn.replace(".gz", ".bgz")
    ungzip(in_fn, tmp_out_fn)
    pysam.tabix_compress(tmp_out_fn, out_fn, force = True)
    os.unlink(tmp_out_fn)
    
    return out_fn

def merge_pileup_tables(pileup_fns, output_pileup_fn, tmp_dir, verbose = True):
    """ sort and merge the sense and antisense pileup tables
    samtools returns pileup format sorted by largest chromosome in 
    descending order
    values for duplicated intervals then can be summed in one pass
    see https://stackoverflow.com/questions/23450145/sort-a-big-file-with-python-heapq-merge
    """
    
    if verbose:
        print("started concatenating anti-sense and sense pileup tabls",
              file = sys.stderr)
              
    ## concat the pileup tables
    out_tmp_ptable = os.path.join(tmp_dir, "tmp_pileup_table.tsv.gz")
    
    nlines = 0
    with gzip.open(out_tmp_ptable, 'wt') as out_tmp_ptable_fn:
      for idx, pfn in enumerate(pileup_fns):
        fn = gzip.open(pfn, 'rt')
        header = fn.readline()
        
        if idx == 0:
          out_tmp_ptable_fn.write(header)
          nlines += 1
        
        for line in fn:
          out_tmp_ptable_fn.write(line)
          nlines += 1

        fn.close()
    
    if verbose:
        print("started sorting anti-sense and sense pileup tables",
              file = sys.stderr)
    
    nlines_per_chunk, n_chunks = compute_chunks(nlines, 1000000)    
    
    if verbose:
        print("chunking pileup tables into {} lines in {} files".format(nlines_per_chunk, n_chunks),
              file = sys.stderr)
    
    ## chunk into sep. files and sort in memory
    chunk_names = []
    with gzip.open(out_tmp_ptable, 'rt') as input_file:
        header = input_file.readline() 
        for chunk_number in itertools.count(1):
            # read in next chunk of lines and sort them
            lines = itertools.islice(input_file, nlines_per_chunk)
            formatted_lines = []
            for x in lines:
                vals = x.split("\t")
                vals[1] = int(vals[1]) # start
                formatted_lines.append(vals)

            sorted_chunk = sorted(formatted_lines, key = itemgetter(0,1,2))
            if not sorted_chunk:
                # end of input
                break
             
            chunk_name = os.path.join(tmp_dir, 'chunk_{}.chk'.format(chunk_number))
            chunk_names.append(chunk_name)
            with gzip.open(chunk_name, 'wt') as chunk_file:
                for line in sorted_chunk:
                  line[1] = str(line[1])
                  chunk_file.write("\t".join(line))
    
    if verbose:
        print("merging sorted tables",
              file = sys.stderr)
              
    with ExitStack() as stack, gzip.open(out_tmp_ptable, 'wt') as output_file:
        files = [stack.enter_context(gzip.open(chunk, 'rt')) for chunk in chunk_names]
        output_file.write(header)
        output_file.writelines(heapq.merge(*files, key = lambda x: (x.split("\t")[0], 
                                                                    x.split("\t")[1],
                                                                    x.split("\t")[2])))
    for i in chunk_names:
      os.unlink(i)
    
    if verbose:
        print("summing up columns from sorted tables",
              file = sys.stderr)
              
    pileup_fn = gzip.open(out_tmp_ptable, 'rt')
    pileup_fout = gzip.open(output_pileup_fn, 'wt')
    header = pileup_fn.readline()
    pileup_fout.write(header)
    
    # group by chrom, pos, and strand and sum values
    pileup_generator = file_to_pileup(pileup_fn)
    for ivl, vals in itertools.groupby(pileup_generator, 
                                       key = lambda x : x.ivl_index()):

        aux_vals = []
        for pileup in vals:
            aux_vals.append(pileup.aux_fields)
            base = pileup.ref_base
        
        if len(aux_vals) > 1:
          summed_vals = [sum(x) for x in zip(*aux_vals)]
        else:
          summed_vals = aux_vals[0]
        
        ivl = list(ivl)
        outline = ivl + [summed_vals[0]] + [base] + summed_vals[1:]
        outline = [str(x) for x in outline]
        pileup_fout.write("\t".join(outline) + "\n")
        
    # clean up files    
    pileup_fn.close()
    pileup_fout.close()
    os.unlink(out_tmp_ptable)
    for f in pileup_fns:
      os.unlink(f)

def format_tbls(fn1, fn2, outfn):
    
    df_pos = pd.read_table(fn1, sep = "\t")
    df_neg = pd.read_table(fn2, sep = "\t")

    df = pd.concat([df_pos, df_neg])
    df = df.drop('start', 1)
    df = df.rename(columns = {'end':'pos'})
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
                        default = 'single')
                        
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
                        default = 5, type = float)
                        
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
                        
    parser.add_argument('-i',
                        '--tabix-index',
                        help=textwrap.dedent("""\
                        If set, report pileup table 
                        as bgzip'd and tabix indexed
                        (default: %(default)s)
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
    
    args = parser.parse_args()
    
    bam_name = args.bam
    
    ### check system exectuables
    
    if not is_tool("samtools"):
        sys.exit("samtools is not in path")
    if not is_tool("mpileupToReadCounts"):
        sys.exit("mpileupToreadcounts is not in path, please add rnastruct/src to your PATH")
    if not is_tool("filterBam"):
        sys.exit("filterBam is not in path, please add rnastruct/src to your PATH")

    # ok to have repeated args in samtools command
    pileup_args =  " --count-orphans -x -d 1000000 -L 1000000 -B " + args.pileup    
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
        
    output_pileup_fn = os.path.join(tmp_dir, "pileup_table.tsv.gz")
    if len(pileup_tbls) > 1:
      # paired end
      # merge_pileups
      merge_pileup_tables(pileup_tbls, output_pileup_fn, tmp_dir, verbose)
      
    elif len(pileup_tbls) == 1:
      # single end
      # just move the file
      shutil.move(pileup_tbls[0], output_pileup_fn)
    
    print("parsing pileup into bedgraph format", file = sys.stderr)

    ## parse output into bedgraphs
    generate_bedgraphs(output_pileup_fn, depth, tmp_dir, threads, nucleotides)
    
    format_tbls(os.path.join(tmp_dir, "pos_pileup_all.tmp"), 
                os.path.join(tmp_dir, "neg_pileup_all.tmp"),
                outpre + "pileup_table.tsv.gz")

    print("merging redundant bedgraph entries", file = sys.stderr)
    
    ## merge redundant intervals
    merge_bedgraphs(tmp_dir, "pos_", outpre)
    merge_bedgraphs(tmp_dir, "neg_", outpre)
    
    ## bgzip and index if requested
    if args.tabix_index:
        out_tbl = outpre + "pileup_table.tsv.gz" 
        out_bgzip = bgzip(out_tbl)        
        pysam.tabix_index(out_bgzip, 
                          seq_col = 0, 
                          start_col = 1,
                          end_col = 1, 
                          zerobased = False,
                          force = True, 
                          line_skip = 1)

                
if __name__ == '__main__': main()

