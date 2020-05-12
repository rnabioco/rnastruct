#!/usr/bin/env python3
import argparse
import textwrap
import os 
import sys
import pysam
import math
import gzip 
import multiprocessing as mp
import pandas as pd
import numpy as np
import shutil
import subprocess
import functools

from extract_transcript import Structpileup
from utils import bgzip
import norm_methods

def unix_sort2(infile, outfile, threads, memory = "8G", header = True, compress = True,
        uncompress = True, verbose = True):
  
    sort_test = subprocess.run(["sort", "--parallel=2"],
              input = "hello world",
              encoding = 'ascii')
    if sort_test.returncode == 0:
      parallel_args = "--parallel=" + str(threads)
    else:
      parallel_args = ""
      if verbose:
        print("Parallel sort not available (man sort), using 1 thread",
               sys.stderr)
    
    outfn = open(outfile, 'w')
    
    if uncompress:
        if header:
            sort_cmd = "(gunzip -c " + infile + " | head -n 1 && gunzip -c " + infile + \
               " | tail -n+2  | sort -S " + memory + parallel_args + " -k1,1 -k2,2n -k3,3 " + \
               " - "
        else:
            sort_cmd = "gunzip -c " + infile + \
               " | sort -S " + memory + parallel_args + " -k1,1 -k2,2n -k3,3 " + \
               " - "
    else:
        if header:
            sort_cmd = "(head -n 1 " + infile + " && tail -n +2 " + infile + \
               " | sort -S " + memory + parallel_args + " -k1,1 -k2,2n -k3,3 " + \
               " - "
        else:
            sort_cmd = "sort -S " + memory + parallel_args + "-k1,1 -k2,2n -k3,3 " + infile 

    if compress and header:
        sort_cmd += " | gzip )"
        sort_cmd = sort_cmd.replace("| head -n 1 &&", "| head -n 1 | gzip &&")
    elif compress:
        sort_cmd += " | gzip"

    if verbose:
        print("sort command is:\n" + sort_cmd, file = sys.stderr)

    sort_run = subprocess.run(sort_cmd, shell=True,
            stderr = sys.stderr,
            stdout = outfn)

    outfn.close()
    return outfile

def parse_tbx_to_pileup(fh, **kwargs):
    for line in fh.fetch(**kwargs):
        row = Structpileup(line.rstrip())
        yield row

def parse_line_to_pileup(fh, header = True):
    
    hdr = fh.readline()
    if header:
        hdr = fh.readline()
    
    for line in fh:
        row = Structpileup(line.rstrip())
        yield row

def combine_std_err(stderr_t, stderr_u):
    return math.sqrt(stderr_t**2 + stderr_u**2)

def np_combine_std_err(stderr_t, stderr_u):
    return math.sqrt(stderr_t**2 + stderr_u**2)

def calc_global_norm(x):
    """ boxplot outlier removal with upper 90% normalization """
    qs = np.quantile(x, [0.25, 0.75])
    iqr = np.diff(qs)[0]
    high_cutoff = qs[1] + (1.5 * iqr)
    
    # based on shapemapper approach, use either outlier or top 10% of
    # data, whichever removes the fewest points

    q = sum(x >= high_cutoff)
    ten_percent = sum(x > np.quantile(x, 0.90))
    
    if (q <= ten_percent):
      low_cutoff = np.quantile(x[x < high_cutoff], 0.90)
      norm_factor = np.mean(x[(x > low_cutoff) & (x < high_cutoff)])
      print("normalization factor = ", 
            norm_factor,
            " using boxplot",
            " with ", q, " removed",
            file = sys.stderr)
    else:
      idx = (x <= np.quantile(x, 0.90)) & ( x > np.quantile(x, 0.80))
      norm_factor = np.mean(x[idx])
      print("normalization factor = ", 
            norm_factor,
            " using top 10% removal ",
            " with ", ten_percent, " removed",
            file = sys.stderr)
    return norm_factor

def remove_untreated(t_path, ut_path, contig, depth, nucs_to_keep,
        outfile, max_untreated_ratio, max_treated_ratio):

    print("processing ", contig)  
    t_df = pd.read_csv(t_path,
            sep = "\t", compression = "gzip") 
    ut_df = pd.read_csv(ut_path,
            sep = "\t", compression = "gzip")
    
    t_df = t_df[t_df.iloc[:, 4] >= depth]

    ut_df = ut_df[['chr', 'pos', 'strand', 'mutation_ratio', 'stderr']]

    out = pd.merge(t_df, ut_df, on = ["chr", "pos", "strand"], how = 'left', 
            suffixes = ("", "_ut"))
    
    out = out.fillna(0.0)
    
    out = out[out.ref_base.isin(nucs_to_keep)]

    # drop ivls with high mutation rate in untreated 
    out = out[out['mutation_ratio_ut'] < max_untreated_ratio]
    
    # drop ivls with high modification rate in treated 
    out = out[out['mutation_ratio'] < max_treated_ratio]
    out = out.reset_index(drop=True)
    
    out["mutation_ratio_bg_sub"] = out["mutation_ratio"] - out["mutation_ratio_ut"]
    
    out["stderr_bg_sub"] = np.sqrt(np.power(out["stderr"], 2) + 
            np.power(out["stderr_ut"], 2))
    
    #out = out.drop(['mutation_ratio_ut', 'stderr_ut'], axis=1)

    out.to_csv(outfile,
            sep = "\t", index = False,  compression = "gzip")
    return outfile


def split_files(contig, t_fn, ut_fn, outdir, header):
    
    t_tbx = pysam.TabixFile(t_fn)
    u_tbx = pysam.TabixFile(ut_fn)
   
    print("splitting ", contig)  
    with gzip.open(os.path.join(outdir, "t_" + contig + ".tsv.gz"), 'wt', compresslevel = 6) as f:
      f.write(header)  
      for ivl in t_tbx.fetch(reference = contig):
         f.write(ivl + "\n")
    
    with gzip.open(os.path.join(outdir, "ut_"+ contig + ".tsv.gz"), 'wt', compresslevel = 6) as f:
      f.write(header)  
      for ivl in u_tbx.fetch(reference = contig):
         f.write(ivl + "\n")
    t_tbx.close()
    u_tbx.close()

    return True

def main():
    
    parser = argparse.ArgumentParser(description="""
    filter and normalize DMS profiles from treated and untreated samples
    """,

    formatter_class = argparse.RawTextHelpFormatter )

    parser.add_argument('-t',
                        '--treated-tabix-file',
                        help ="""tabix indexed pileup file
                        \n""",
                        required = True)

    parser.add_argument('-u',
                        '--untreated-tabix-file',
                        help ="""tabix indexed pileup file
                        \n""",
                        required = True)
    
    parser.add_argument('-o',
                        '--outpre',
                        help ="""tabix indexed outfile prefix
                        \n""",
                        required = True)

    parser.add_argument('-ur',
                        '--untreated_reactivity_cutoff',
                        help = textwrap.dedent("""\
                        reactivities cutoffs to exclude.
                        (default: %(default)s)
                        \n"""),
                        required = False,
                        default = 0.02,
                        type = float)
    
    parser.add_argument('-tr',
                        '--treated_reactivity_cutoff',
                        help = textwrap.dedent("""\
                        reactivities cutoffs to exclude.
                        (default: %(default)s)
                        \n"""),
                        required = False,
                        default = 0.10,
                        type = float)
    
    parser.add_argument('-d',
                        '--depth_cutoff',
                        help = textwrap.dedent("""\
                        depth cutoffs to exclude.
                        (default: %(default)s)
                        \n"""),
                        required = False,
                        default = 10,
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
                        default = "ATCG",
                        type = str)
    
    parser.add_argument('-p',
                        '--proc',
                        help = textwrap.dedent("""\
                        n processors
                        \n"""),
                        required = False,
                        default = 1,
                        type = int)

    args = parser.parse_args() 
    
    verbose = True
    outdir = "."
    thread_count = args.proc
    memory_chunk = "1G"
    depth = args.depth_cutoff

    with gzip.open(args.treated_tabix_file, 'rt') as f:
      hdr = f.readline()
      f.seek(0)
    
    outpre = args.outpre
    outdir = os.path.dirname(outpre)
    if not os.path.exists(outdir) and outdir:
       os.makedirs(outdir)

    tmp_dir = os.path.join(outdir, "tmpfiles")
    if not os.path.exists(tmp_dir):
       os.makedirs(tmp_dir)
       if verbose:
           print("temporary files placed in directory:\n{}".format(tmp_dir),
                 file = sys.stderr)
       else:
         sys.exit("temporary files directory already exists:\n{}".format(tmp_dir))

  
    t_tbx = pysam.TabixFile(args.treated_tabix_file)
    u_tbx = pysam.TabixFile(args.untreated_tabix_file)
    contigs = t_tbx.contigs
    
    # split into per chrom files
    
    t_tbx_fn = args.treated_tabix_file
    ut_tbx_fn = args.untreated_tabix_file
    
    
    func = functools.partial(split_files,
                             t_fn = t_tbx_fn,
                             ut_fn = ut_tbx_fn,
                             outdir = tmp_dir,
                             header = hdr)


    with mp.Pool(thread_count) as p:
        chunked_files = p.map(func, contigs)

    # operate over each contig df in parallel
    t_contig_files = [os.path.join(tmp_dir, "t_" + contig + ".tsv.gz") for
            contig in contigs]
    ut_contig_files = [os.path.join(tmp_dir, "ut_" + contig + ".tsv.gz")
            for contig in contigs]
    out_contig_files = [os.path.join(tmp_dir, "merge_" + contig +
            ".tsv.gz") for contig in contigs]
    
    nucs_to_keep = args.nucleotides
    nucs_to_keep = [x.upper() for x in nucs_to_keep]

    mp_args = zip(t_contig_files, 
                  ut_contig_files, 
                  contigs, 
                  [depth] * len(contigs), 
                  [nucs_to_keep] * len(contigs), 
                  out_contig_files,
                  [args.untreated_reactivity_cutoff] * len(contigs),
                  [args.treated_reactivity_cutoff] * len(contigs))

    with mp.Pool(thread_count) as pool:
        contig_files = pool.starmap(remove_untreated, mp_args) 

    output = os.path.join(tmp_dir, "pileup_table.tsv.gz")

    # combine gzipped per chromosome files
    output_tbl = gzip.open(output, 'wt', compresslevel = 6)
    for idx, fn in enumerate(contig_files):
        with gzip.open(fn, 'rt') as data:

            if idx == 0:
                # keep header from first file
                shutil.copyfileobj(data, output_tbl)
            else:
                # skip header for remaining
                data.readline()
                shutil.copyfileobj(data, output_tbl)
    output_tbl.close()
    
    # calc norm factor.
    df = pd.read_csv(output, sep = "\t", compression = 'gzip')
    #normalization_factor = calc_global_norm(df["mutation_ratio_bg_sub"].values)
    
    normalization_factor = np.quantile(df["mutation_ratio_bg_sub"].values, [0.95])
    windsorized_vals = norm_methods.windsorize(df["mutation_ratio_bg_sub"].values)

    df["mutation_ratio_norm"] = windsorized_vals / normalization_factor 
    df["stderr_norm"] = df["stderr_bg_sub"] / normalization_factor
    df.to_csv(output, index = False, sep = "\t", compression = "gzip")
    
    print(outdir)
    # sort output
    sorted_output = outpre + "pileup_table.tsv.gz"
    unix_sort2(output, sorted_output, memory = memory_chunk, threads = thread_count)
    
    # bgzip output
    out_bgzip = bgzip(sorted_output)
    pysam.tabix_index(out_bgzip,
                      seq_col = 0,
                      start_col = 1,
                      end_col = 1,
                      zerobased = False,
                      force = True,
                      line_skip = 1)

    # cleanup
    debug = False
    if not debug:
        os.unlink(output)
        
        tmpfiles = contig_files + t_contig_files + ut_contig_files
        for fn in tmpfiles:
            os.unlink(fn)
        os.rmdir(tmp_dir)

if __name__ == '__main__': main()


