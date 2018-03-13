#!/usr/bin/env python3
import argparse
import textwrap
import os 
import sys
import subprocess
import gzip
import pandas as pd
import uuid 
import shutil
import itertools
import heapq
import io
import multiprocessing as mp
from contextlib import ExitStack
from datetime import datetime
from collections import Counter
from functools import partial
from operator import itemgetter

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

def line_to_interval(line):
    line = line.rstrip()
    fields = line.split("\t")
    return Interval(fields[0], 
                    int(fields[1]), 
                    int(fields[2]),
                    float(fields[3]))

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
            outbg.write(str(ivl_cache) + "\n")
            ivl_cache = ivl
            current_chrom = ivl.chrom
        elif ivl_cache.end == ivl.start and ivl_cache.count == ivl.count:
            # same coverage, modify end
            ivl_cache.end = ivl.end
        else:
            outbg.write(str(ivl_cache) + "\n")
            ivl_cache = ivl
            
    ## clear out the last interval        
    outbg.write(str(ivl_cache) + "\n")    
    
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

def merge_bedgraphs(prefix, strand, 
                    insuffix = ".bedgraph.tmp",
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
    
def format_bedgraphs(df, depth, prefix):
    
    """ take pandas dataframe and split into bedgraph dataframs
    return list of dataframes and list of output filenames
    pass to write bedgraphs.
    
    Note use of global variable NUCS here to denote if bedgraphs should be
    restricted to A and C bases only
    """
    df = df.assign(mismatch_ratio = lambda df: df.mmcount / df.depth)
    df = df.assign(insertion_ratio = lambda df: df.inscount / df.depth)
    df = df.assign(deletion_ratio = lambda df: df.delcount / df.depth)

    df = df.assign(start = lambda df:df.pos - 1)
    df = df.rename(columns = {'pos':'end'})

    df = df[df.depth >= depth]
    df_depth = df[['chr', 'start', 'end', 'depth']]
    
    df = df[df.ref_base.isin(NUCS)]
    
    df_mm = df[['chr', 'start', 'end', 'mismatch_ratio']]
    df_ins = df[['chr', 'start', 'end', 'insertion_ratio']]
    df_del = df[['chr', 'start', 'end', 'deletion_ratio']]

    df_mm = df_mm.sort_values(['chr', 'start'], ascending=[True, True])
    df_ins = df_ins.sort_values(['chr', 'start'], ascending=[True, True])
    df_del = df_del.sort_values(['chr', 'start'], ascending=[True, True])
    df_depth = df_depth.sort_values(['chr', 'start'], ascending = [True, True])
    
    ## run merge intervals on memory mapped data
    df_list = [df_mm, df_ins, df_del, df_depth]
    df_fnames = [
      "mismatches",
      "insertions",
      "deletions",
      "depth"
    ]
    
    mem_mapped_dfs = []
    mem_mapped_dfs_fnames = []
    for idx, df in enumerate(df_list):
        if df.empty:
          continue
        
        ## convert to in-memory file and merge intervals in place
        fobj = memmap_df(df)
        mem_mapped_dfs.append(fobj)
        mem_mapped_dfs_fnames.append(df_fnames[idx])
    
    out_fns = [prefix + x + ".bedgraph.tmp" for x in mem_mapped_dfs_fnames]
    
    return(mem_mapped_dfs, out_fns)

def parse_library_type(bam_path, strandedness, libtype):
    """ return a list of appropriate flags for filtering bamfile 
      returns a list with a list for each bam flag setting
      first element is a filtering flag
      second element is a string that indicates if the alignments are sense or antisense
    """
    
    if strandedness == "fr-firststrand":
      if libtype == "paired":
        ## parse out sense alignments (R2) (second in pair)
        sam_flag_1 = "-f 128 "
        lib_type_1 = "sense"
        ## parse out antisense alignments (R1) (first in pair)
        sam_flag_2 = "-f 64 "
        lib_type_2 = "antisense"
        
        return [[sam_flag_1, lib_type_1], [sam_flag_2, lib_type_2]]
        
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
        
        return [[sam_flag_1, lib_type_1], [sam_flag_2, lib_type_2]]
      else:
        # single end
        return [["", "sense"]] 
        
    elif strandedness == "unstranded":
      
      return [["", "sense"]]
      
    else:
      sys.exit("libtype specification failed")
    
def generate_mismatch_profile(input_bam, fasta, additional_args, depth, outpre, 
        threads, deletion_length, bam_flag, libtype, debug = False):
    
   if debug:
       print("started processing {}".format(str(datetime.now())),
           file = sys.stderr)

   ## generate per nucleotide mismatch and indel counts
   if threads == 1:
       output = generate_pileup(input_bam, fasta, depth, 
                                deletion_length, additional_args,
                                bam_flag, libtype,
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
       new_pres = [os.path.join(outpre, x + "_") for x in contigs] 
       
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
   
   return(output)

def split_and_apply(df, min_depth, outprefix):
    """ parse into strand-specific bedgraphs (pos and neg)
    and return a list of bedgraph dataframes and their output filenames.
    pass to write_bedgraph to write to disk
    """
    
    df_pos = df[df.strand == "+"]
    df_neg = df[df.strand == "-"]
    
    bgs_pos, bgs_pos_fns = format_bedgraphs(df_pos, min_depth, outprefix + "pos_") 
    bgs_neg, bgs_neg_fns = format_bedgraphs(df_neg, min_depth, outprefix + "neg_")
    
    bgs_out = bgs_pos + bgs_neg
    bgs_fns = bgs_pos_fns + bgs_neg_fns
    
    return([bgs_out, bgs_fns])
    
def generate_bedgraphs(pileup_fn, depth, outpre, threads, chunk_size = 100000):
   """ master function for generating bedgraphs in parallel """
   
   ## read in pileup table in chunks to keep memory low
   reader = pd.read_table(pileup_fn, compression = 'gzip', chunksize = chunk_size)
   
   ## apply multithreading if cpus available
   pool = mp.Pool(threads)
   
   func = partial(split_and_apply,
            min_depth = depth,
            outprefix = outpre)
  
   ## results is a list of pandas df's and output filenames
   results = pool.imap(func, reader)
   pool.close()
   pool.join()
   
   for res in results:
       write_bedgraphs(res[0], res[1])
   
   ## merge redundant intervals that may flank each chunk
   merge_bedgraphs(outpre, "pos_")
   merge_bedgraphs(outpre, "neg_")
   
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
    out_tmp_ptable = os.path.join(tmp_dir, "_tmp_pileup_table.tsv.gz")
    with gzip.open(out_tmp_ptable, 'wt') as out_tmp_ptable_fn:
      for idx, pfn in enumerate(pileup_fns):
        fn = gzip.open(pfn, 'rt')
        header = fn.readline()
        if idx == 0:
          out_tmp_ptable_fn.write(header)
        for line in fn:
          out_tmp_ptable_fn.write(line)
        fn.close()
    
    if verbose:
        print("started sorting anti-sense and sense pileup tables",
              file = sys.stderr)
    
    ## chunk into sep. files and sort in memory
    chunk_names = []
    with gzip.open(out_tmp_ptable, 'rt') as input_file:
        header = input_file.readline() 
        for chunk_number in itertools.count(1):
            # read in next 100k lines and sort them
            lines = itertools.islice(input_file, 100000)
            lines = [x.split("\t") for x in lines]
            sorted_chunk = sorted(lines, key = itemgetter(0,1,2))
            if not sorted_chunk:
                # end of input
                break
    
            chunk_name = os.path.join(tmp_dir, 'chunk_{}.chk'.format(chunk_number))
            chunk_names.append(chunk_name)
            with gzip.open(chunk_name, 'wt') as chunk_file:
                for line in sorted_chunk:
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
                        mismatc and indel ratios. Provide
                        as a string. i.e. to report 
                        for all nucleotides. "ATCG"
                        (default: %(default)s)
                        \n"""),
                        required = False,
                        default = "AC",
                        type = str)
                        
    parser.add_argument('-v',
                        '--verbose',
                        help="""print run information (default: %(default)s)\n""",
                        type = bool,
                        required = False,
                        default = False)

    args = parser.parse_args()
    
    bam_name = args.bam

    # ok to have repeated args in samtools command
    pileup_args =  " --count-orphans -x -d 1000000 -L 1000000 -B " + args.pileup    
    fasta_name = args.fasta
    depth = args.depth
    outpre = args.outpre
    deletion_length = args.deletion_length
    threads = args.threads
    verbose = args.verbose
    library = args.library
    nucs = args.nucleotides
    strandedness = args.strandedness
    
    #### set up global
    global NUCS
    NUCS = list(nucs)
    
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

    #### parse library type  
    bam_flags = parse_library_type(bam_name, strandedness, library)
    
    ## parse bam into two new bams if paired end
    ## one bam with all alignments that report the correct strand of the fragment
    ## the other bam with all alignments that are rev-comp of the fragment
    ## invert the second bam reported strands and merge
    pileup_tbls = []
      
    for bam_flag in bam_flags:
      bam_flag_filter = bam_flag[0]
      align_type = bam_flag[1]
      pileup_fn = os.path.join(tmp_dir, bam_flag[1])
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
        verbose) 
      
      pileup_tbls.append(pileup_tbl)
        
    output_pileup_fn = outpre + "pileup_table.tsv.gz"
    if len(pileup_tbls) == 2:
      # paired end
      # merge_pileupes
      merge_pileup_tables(pileup_tbls, output_pileup_fn, tmp_dir, verbose)
      
    elif len(pileup_tbls) == 1:
      # single end
      # just move the file
      shutil.move(pileup_tbls[0], output_pileup_fn)
   
    ## parse output into bedgraphs
    generate_bedgraphs(output_pileup_fn, depth, outpre, threads)
    
    print("removing temp directory: {}".format(tmp_dir))
    shutil.rmtree(tmp_dir, ignore_errors=False, onerror=None)     

                
if __name__ == '__main__': main()

