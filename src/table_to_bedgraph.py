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
from count_mutations import Interval

class StrandedInterval:
    """ stranded bed interval object """
    def __init__(self, chrom, start, end, strand, values):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.vals = values
    def __str__(self):
        out_vals = [str(x) for x in self.vals]
        return "{}\t{}\t{}\t{}".format(self.chrom,
                self.start,
                self.end,
                "\t".join(out_vals))


def tbxline_to_interval(line, pos_cols, cols, sc_col):

    line = line.rstrip()
    fields = line.split("\t")
    aux = [fields[i] for i in cols]
    
    if pos_cols[0] == pos_cols[1]:
      end = int(fields[pos_cols[0]])
      start = end - 1
    else:
      start = int(fields[pos_cols[0]])
      end = int(fields[pos_cols[1]])
    
    return StrandedInterval(fields[0], 
                    start, 
                    end,
                    fields[sc_col],
                    aux)

        
def convert_tabix_to_bg(bg, pos_cols, cols, strand, sc_col, outbg, header = True):
    """ convert pileup format to bedgraph 
    merging adjacent intervals with the same
    count values
    """
    if header:
        hdr = bg.readline()
    
    strand_incorrect = True
    while strand_incorrect:
        ivl_cache = tbxline_to_interval(bg.readline(), pos_cols, cols,
                sc_col)
        current_chrom = ivl_cache.chrom
        if ivl_cache.strand == strand:
            strand_incorrect = False

    for line in bg:
        ivl = tbxline_to_interval(line, pos_cols, cols, sc_col)
        if ivl.strand != strand:
            continue

        if ivl.chrom != current_chrom:
            outbg.write("{}\n".format(str(ivl_cache)))
            ivl_cache = ivl
            current_chrom = ivl.chrom
        elif ivl_cache.end == ivl.start and ivl_cache.vals == ivl.vals:
            # same value and adjacent, modify end
            ivl_cache.end = ivl.end
        else:
            outbg.write("{}\n".format(str(ivl_cache)))
            ivl_cache = ivl
            
    ## clear out the last interval        
    outbg.write("{}\n".format(str(ivl_cache)))
    

def main():
    
    parser = argparse.ArgumentParser(description="""
    extract regions from a tabix file and convert to a bedgraph
    """,
    formatter_class = argparse.RawTextHelpFormatter )

    parser.add_argument('-i',
                        '--input_table',
                        help ="""tabix indexed input
                        \n""",
                        required = True)
    parser.add_argument('-b',
                        '--beginning_col',
                        help = textwrap.dedent("""\
                        start column (one based)
                        \n"""), 
                        required = True)
    parser.add_argument('-e',
                        '--end_col',
                        help = textwrap.dedent("""\
                        end column (one based), 
                        if end column is same as start column, 
                        will be interpeted as one-based position
                        \n"""), 
                        required = True)
    parser.add_argument('-sc',
                        '--strand_col',
                        help = textwrap.dedent("""\
                        column with strand info (one-based)
                        \n"""), 
                        required = True)
    parser.add_argument('-c',
                        '--cols',
                        help = textwrap.dedent("""\
                        cols to keep 
                        \n"""), 
                        required = False, 
                        default = [24, 25], 
                        nargs = "+")
    
    parser.add_argument('-s',
                        '--strand',
                        help = textwrap.dedent("""\
                        strand to return
                        \n"""), 
                        required = True, 
                        default = "+")
    
    parser.add_argument('-n',
                        '--noheader',
                        help = textwrap.dedent("""\
                        if set, indicats that input pileup table does not 
                        have a header
                        \n"""), 
                        required = False, 
                        action = 'store_false')

    args = parser.parse_args()
    
    in_fn = args.input_table
    
    s_col = int(args.beginning_col) - 1
    e_col = int(args.end_col) - 1
    sc_col = int(args.strand_col) - 1
    cols = [int(x) - 1 for x in args.cols] 
    in_fh = gzip.open(in_fn, 'rt')
    
    convert_tabix_to_bg(in_fh, [s_col, e_col], cols, args.strand, sc_col, sys.stdout, header = args.noheader)

                
if __name__ == '__main__': main()

