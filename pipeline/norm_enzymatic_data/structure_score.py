#! /usr/bin/env python3

import argparse
import gzip as gz
from math import *
import tempfile 
import sys
import binascii

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
    return Interval(fields[0], 
                    int(fields[1]), 
                    int(fields[2]),
                    [float(x) for x in fields[3:]])

def merge_ivls(current_ivl, previous_ivl):
  previous_ivl.end = current_ivl.end
  return previous_ivl

def ivl_is_same(current_ivl, previous_ivl):
  if current_ivl.chrom != previous_ivl.chrom:
      return False
  if current_ivl.start != previous_ivl.end:
      return False
  if current_ivl.vals != previous_ivl.vals:
      return False
  else:
      return True

def is_gz_file(filepath):
    """
    "https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed"
    """
    with open(filepath, 'rb') as test_f:
        return binascii.hexlify(test_f.read(2)) == b'1f8b'

def gzopen(fn, *args):
    
    if(is_gz_file(fn)):
      f = gz.open(fn, *args)
    else:
      f = open(fn, *args)
    
    return f

def count_nascent_norm(nascent_count,
               ss_count,
               ds_count,
               read_counts):
    """
    returns tuple with normalized nascent, normalized_ss,
    and normalized ds.

    Tassa's method
    """

    norm_nascent = 1e6 * (nascent_count / float(read_counts[0]))
    if norm_nascent == 0:
        return 0, None, None

    norm_ss = (1e6 * (ss_count / float(read_counts[1]))) / norm_nascent 
    norm_ds = (1e6 * (ds_count / float(read_counts[2]))) / norm_nascent
    
    return norm_nascent, norm_ss, norm_ds

def count_norm(nascent_count,
               ss_count,
               ds_count,
               read_counts):
    """
    returns tuple with normalized nascent, normalized_ss,
    and normalized ds
    the first element of read_counts will be used to normalize
    nascent RNA (for now)
    approach from here:
    https://dx.doi.org/10.1105%2Ftpc.112.104232
    """

    norm_nascent = 1e6 * (nascent_count / float(read_counts[0]))
    
    if ss_count == 0 and ds_count == 0:
        return norm_nascent, None, None
    
    norm_ss = ss_count * \
            (max(read_counts[1], read_counts[2]) /
             float(read_counts[1]))

    norm_ds = ds_count * \
            (max(read_counts[1], read_counts[2]) /
             float(read_counts[2]))
    
    return norm_nascent, norm_ss, norm_ds

def coverage_norm(nascent_count,
               ss_count,
               ds_count,
               cov_lengths):
    """
    returns tuple with normalized nascent, normalized_ss,
    and normalized ds
    the first element of cov_lengths will be used to normalize
    nascent RNA (for now)
    approach from here:
    https://doi.org/10.1016/j.molcel.2014.12.004
    """

    norm_nascent = 1e6 * (nascent_count / float(cov_lengths[0]))
    
    if ss_count == 0 and ds_count == 0:
        return norm_nascent, None, None
    
    norm_ss = ss_count * \
            (max(cov_lengths[1], cov_lengths[2]) /
             float(cov_lengths[1]))

    norm_ds = ds_count * \
            (max(cov_lengths[1], cov_lengths[2]) /
             float(cov_lengths[2]))
    
    return norm_nascent, norm_ss, norm_ds

def calc_score(in_fn, read_cnts, score_fields, norm_method):
  norm_methods = ["counts", "coverage", "nascent"]

  if norm_method not in norm_methods:
      sys.exit("unknown normalization method: {}".format(norm_method))

  with gzopen(in_fn,'rt') as fin:        
    previous_ivl = None
    
    for idx,line in enumerate(fin):        
        if idx == 0:
          previous_ivl = line_to_interval(line)
          ivl = previous_ivl
        else:
          ivl = line_to_interval(line)
        
        nascent_score = ivl.vals[score_fields[0]]
        ss_score = ivl.vals[score_fields[1]]
        ds_score = ivl.vals[score_fields[2]]

        if norm_method == "counts":
            norm_nascent, norm_ss, norm_ds = count_norm( nascent_score, ss_score, ds_score, read_cnts)

        elif norm_method == "coverage":
            norm_nascent, norm_ss, norm_ds = coverage_norm( nascent_score, ss_score, ds_score, read_cnts)
        
        elif norm_method == "nascent":
            norm_nascent, norm_ss, norm_ds = count_nascent_norm( nascent_score, ss_score, ds_score, read_cnts)

        # don't report intervals if normalization results in NaN or if
        # there is no ss or ds seq coverage 
        if norm_ss is None:
            if previous_ivl is not None and idx != 0:
                print(previous_ivl)
            previous_ivl = None
            continue
        
       
        glog_norm_struct = \
                log(norm_ds + sqrt(1 + norm_ds ** 2), 2) - \
                log(norm_ss + sqrt(1 + norm_ss ** 2), 2)
        
        ivl.vals += [norm_nascent, norm_ss, norm_ds, glog_norm_struct]
        fields = ivl.vals
        fields = [str(x) for x in fields] 
        
        if previous_ivl is None:
            previous_ivl = ivl
        elif ivl_is_same(ivl, previous_ivl):
            previous_ivl = merge_ivls(ivl, previous_ivl)
        elif idx == 0:
            previous_ivl.vals = ivl.vals
        else:
            print(previous_ivl)
            previous_ivl = ivl

    if previous_ivl is None:
        pass
    elif ivl_is_same(ivl, previous_ivl):
        previous_ivl = merge_ivls(ivl, previous_ivl)
        print(previous_ivl)
    else:
        print(previous_ivl)

def main():
  
    parser = argparse.ArgumentParser(description="""
    Calculate a structure score for PIP-seq dataset""",
    formatter_class = argparse.RawTextHelpFormatter )

    parser.add_argument('-c',
                        '--counts',
                        help ="""counts bed-like file""",
                        required = True)

    parser.add_argument('-r',
                        '--read_counts',
                        help ="""csv list of number of reads to normalize
                        against for each library. 
                        first field is total rna,
                        second field is ssRNA,
                        third field is dsRNA.
                        default = 1e6,1e6,1e6
                        """,
                        required = False,
                        default = "1e6,1e6,1e6")

    parser.add_argument('-f',
                        '--fields',
                        help = """csv list of fields to process
                        first field is total rna,
                        second field is ssRNA,
                        third field is dsRNA. 
                        Default = 4,5,6""",
                        required = False,
                        default = "4,5,6")

    parser.add_argument('-m',
                        '--method',
                        help = """
                        Normalization method. 
                        one of "nascent", "counts" or "coverage"
                        """,
                        required = False,
                        default = "nascent")

    args = parser.parse_args()
    

    fields = [int(x) - 4 for x in args.fields.split(",")]
    read_counts = [int(float(x)) for x in args.read_counts.split(",")]
    norm_method = args.method

    calc_score(args.counts, read_counts, fields, norm_method)
if __name__ == '__main__': main()

