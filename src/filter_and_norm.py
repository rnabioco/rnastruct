#!/usr/bin/env python3
import argparse
import textwrap
import os 
import sys
import pysam
import math
import gzip 
import multiprocessing as mp
import uuid
import pandas as pd
import numpy as np
import shutil
import subprocess

from extract_transcript import Structpileup
from count_mutations import bgzip

def unix_sort2(infile, outfile, threads, memory = "8G", header = True, compress = True,
        uncompress = True, verbose = True):
    
    outfn = open(outfile, 'w')
    
    if uncompress:
        if header:
            sort_cmd = "(gunzip -c " + infile + " | head -n 1 && gunzip -c " + infile + \
               " | tail -n+2  | sort -S " + memory + " --parallel=" + str(threads) + " -k1,1 -k2,2n -k3,3 " + \
               " - "
        else:
            sort_cmd = "gunzip -c " + infile + \
               " | sort -S " + memory + " --parallel=" + str(threads) + " -k1,1 -k2,2n -k3,3 " + \
               " - "
    else:
        if header:
            sort_cmd = "(head -n 1 " + infile + " && tail -n +2 " + infile + \
               " | sort -S " + memory + " --parallel=" + str(threads) + " -k1,1 -k2,2n -k3,3 " + \
               " - "
        else:
            sort_cmd = "sort -S " + memory + " --parallel=" + str(threads) + "-k1,1 -k2,2n -k3,3 " + infile 

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

    q = sum(x < high_cutoff)
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

def remove_untreated(t_path, ut_path, contig, depth, outfile):

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

    # drop ivls with high mutation rate in untreated 
    out = out[out['mutation_ratio_ut'] < 0.02]
    out = out.reset_index(drop=True)
    
    out["mutation_ratio_bg_sub"] = out["mutation_ratio"] - out["mutation_ratio_ut"]
    
    out["stderr_bg_sub"] = np.sqrt(np.power(out["stderr"], 2) + 
            np.power(out["stderr_ut"], 2))
    
    #out = out.drop(['mutation_ratio_ut', 'stderr_ut'], axis=1)

    out.to_csv(outfile,
            sep = "\t", index = False,  compression = "gzip")
    return outfile

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
                        '--outdir',
                        help ="""tabix indexed outfile directory
                        \n""",
                        required = True)

    parser.add_argument('-r',
                        '--reactivity_cutoff',
                        help = textwrap.dedent("""\
                        reactivities cutoffs to exclude.
                        (default: %(default)s)
                        \n"""),
                        required = False,
                        default = 0.02,
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
 
    #tmp_dir = os.path.join(outdir, "tmpfiles-" + str(uuid.uuid4()))
    tmp_dir = args.outdir
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
    
    #for contig in contigs:
    #  print("processing ", contig)  
    #  with gzip.open(os.path.join(tmp_dir, "t_" + contig + ".tsv.gz"), 'wt', compresslevel = 6) as f:
    #    f.write(hdr)  
    #    for ivl in t_tbx.fetch(reference = contig):
    #       f.write(ivl + "\n")
    #  
    #  with gzip.open(os.path.join(tmp_dir, "ut_"+ contig + ".tsv.gz"), 'wt', compresslevel = 6) as f:
    #    f.write(hdr)  
    #    for ivl in u_tbx.fetch(reference = contig):
    #       f.write(ivl + "\n")
    
    # operate in parallel
    t_contig_files = [os.path.join(tmp_dir, "t_" + contig + ".tsv.gz") for
            contig in contigs]
    ut_contig_files = [os.path.join(tmp_dir, "ut_" + contig + ".tsv.gz")
            for contig in contigs]
    out_contig_files = [os.path.join(tmp_dir, "merge_" + contig +
            ".tsv.gz") for contig in contigs]

    mp_args = zip(t_contig_files, ut_contig_files, contigs, [depth] * len(contigs), out_contig_files)

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
    normalization_factor = calc_global_norm(df["mutation_ratio_bg_sub"].values)
    df["mutation_ratio_norm"] = df["mutation_ratio_bg_sub"] / normalization_factor 
    df["stderr_norm"] = df["stderr_bg_sub"] / normalization_factor
    df.to_csv(output, index = False, sep = "\t", compression = "gzip")
    
    # sort output
    sorted_output = os.path.join(tmp_dir, "pileup_table.sorted.tsv.gz")
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

    os.unlink(output)
 
if __name__ == '__main__': main()
