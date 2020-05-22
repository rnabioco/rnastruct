#! /usr/bin/env python3
import pysam
import numpy as np
import argparse
import sys
import os
import sys
import gzip
import math
from collections import Counter

import utils
import pdb
import re

bp_comp = {
        "A" : "T",
        "T" : "A",
        "G" : "C",
        "C" : "G",
        "N" : "N"
        }

ref_nts = ["A", "T", "C", "G"]

def allele_counter():
  
    alleles_to_count = ["A","T","C","G"]
    alleles_dict = Counter()
    
    for base in alleles_to_count:
        alleles_dict[base] = 0  
    
    return alleles_dict
            
class IndelCache():
    """cache to hold records for indels. Used to ensure
    that indel position reporting is left aligned if forward stranded
    alignments were used and right aligned if reverse stranded.
    Thinking that RT will hit an adduct from 3'-> 5'. For ambigious
    indels, there isn't a way to precisley localize, however
    left alignment is used in some publications"""

    def __init__(self, chrom = None):
        self.d = Counter()
        self.refdepth = Counter()
        self.chrom = chrom
        
    def add(self, plpc):
        
        for k,v in plpc.indeldict.items():
          self.d[k] += v
        
        for k,v, in plpc.indelcov.items():
          self.refdepth[k] += v
              

    def n(self):
        return len(self.d), len(self.refdepth)

    def get_indel_depth(self, plpc):
        return self.d.pop(plpc.pos, 0)
    
    def get_indel_skipped_coverage(self, plpc):
        return self.refdepth.pop(plpc.pos, 0)
    
    def clear(self):
        self.d.clear()
        self.refdepth.clear()
        self.chrom = None

    def __str__(self):
        res = []
        for k,v in self.d.items():
          out = "{}\t{}".format(str(k), str(v))
          res.append(out)
        return "\n".join(res)


class PileupCounter():
    "Single sample allele counter for samtools mpileup"
    
    def __init__(self, line, strand = "+", max_del = 4):
        self.pos = None
        self.chrom = None
        self.strand = strand
        self.max_del = max_del
        self.base = None
        self.depth = None
        self.plp = None
        self.refcount = 0
        self.mmcount = 0
        self.ncount = 0
        self.mismatch_ratio = 0.0
        self.indel_ratio = 0.0
        self.mutation_ratio = 0.0
        self.stderr = 0.0
        self.indelcount = 0
        self.indelcovcount = 0
        self.indeldict = Counter()
        self.indelcov = Counter()
        self.ac = allele_counter()
        self.read_line(line)
                
    def __str__(self):
        return "\t".join([str(x) for x in self.to_counts()]) + "\n"

    def get_depth(self):
        return self.refcount + self.mmcount + self.indelcovcount

    def mismatch_stats(self):
        depth = self.get_depth() 
        self.mismatch_ratio = float(self.mmcount) / depth
        self.indel_ratio = float(self.indelcount) / depth
        self.mutation_ratio = float((self.mmcount + self.indelcount)) / depth
        self.stderr = math.sqrt(self.mutation_ratio) / math.sqrt(depth)
    
    def to_counts(self):
        indel_counts = self.indelcount
        ref_counts = self.refcount
        mm_counts = self.mmcount
        depth_out = self.get_depth() 
        self.mismatch_stats()
        if self.strand == "-":
            ref_nt = bp_comp[self.base]

            out = [self.chrom,
            self.pos,
            self.strand,
            ref_nt,
            depth_out,
            ref_counts,
            self.ac[bp_comp["A"]],
            self.ac[bp_comp["C"]],
            self.ac[bp_comp["G"]],
            self.ac[bp_comp["T"]],
            mm_counts,
            indel_counts,
            self.mismatch_ratio,
            self.indel_ratio,
            self.mutation_ratio,
            self.stderr
            ]
        else:
            out = [self.chrom,
            self.pos,
            self.strand,
            self.base,
            depth_out,
            ref_counts,
            self.ac["A"],
            self.ac["C"],
            self.ac["G"],
            self.ac["T"],
            mm_counts,
            indel_counts,
            self.mismatch_ratio,
            self.indel_ratio,
            self.mutation_ratio,
            self.stderr
            ]
        return out

    def read_line(self, line):
        fields = line.rstrip().split("\t")
        self.chrom = fields[0]
        self.pos = int(fields[1])
        self.base = fields[2]
        self.depth = int(fields[3])
        self.plp = fields[4]
        self.parse_pileup()
        self.set_depths()
        
    def set_depths(self):
        self.mmcount += sum(self.ac.values())
        self.ac[self.base] += self.refcount
        self.has_indel = len(self.indeldict) > 0

    def parse_pileup(self):
        if self.plp is not None:
            x = 0
            while x < len(self.plp): 
                indelsize_string = ""
                i = self.plp[x]

                if i == "." or i == ",":
                    self.refcount += 1
                    x += 1
                    continue
                
                elif i == "a" or i == "A":
                    self.ac["A"] += 1
                    x += 1
                    continue
                
                elif i == "t" or i == "T":
                    self.ac["T"] += 1
                    x += 1
                    continue
                
                elif i == "c" or i == "C":
                    self.ac["C"] += 1
                    x += 1
                    continue
    
                elif i == "g" or i == "G":
                    self.ac["G"] += 1
                    x += 1
                    continue
    
                elif i == "n" or i == "N":
                    self.ncount += 1
                    x += 1
                    continue
    
                elif i == "<" or i == ">":
                    x += 1
                    continue
    
                elif i == "*":
                    x += 1
                    continue
    
                elif i == "+":
                    # move past + 
                    x += 1
                    while ord(self.plp[x]) >= 48 and ord(self.plp[x]) <= 57 :
                        indelsize_string += self.plp[x]
                        x += 1
                    indel_len = int(indelsize_string)

                    if self.strand == "+":
                        self.indeldict[self.pos] += 1
                    else:
                        self.indeldict[self.pos + 1] += 1

                    x += indel_len
                    continue

                elif i == "-":
                    # move past - 
                    x += 1
                    while ord(self.plp[x]) >= 48 and ord(self.plp[x]) <= 57 :
                        indelsize_string += self.plp[x]
                        x += 1
                    indel_len = int(indelsize_string)
                    if indel_len >= self.max_del:
                      x += indel_len
                      continue

                    if self.strand == "+":
                        self.indeldict[self.pos + indel_len] += 1
                    else:
                        self.indeldict[self.pos + 1] += 1
                    
                    for i in range(self.pos + 1, self.pos + 1 + indel_len):
                      self.indelcov[i] += 1
    
                    x += indel_len
                    continue

                elif i == "$":
                    x += 1
                    continue
                elif i == "^":
                    x += 2 # skip quality
                    continue
                else:
                    sys.exit("unknown char in pileup: {}:{}\n{}".format(x,
                        i, 
                        self.plp))

def write_header(fo):
  
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
  
  fo.write("\t".join(tbl_cols) + "\n")

def mpileup_to_counts(fh, out_fh, min_depth = 10, return_comp = False,
        debug = True, max_del = 4, region = None, header = False): 
      
  if return_comp:
      strand = "-"
  else:
      strand = "+"
  
  if header:
    write_header(out_fh)
  n_good_variants = 0

  ic = IndelCache()
  previous_rec = None
  for i,v in enumerate(fh):
      
      current_rec = PileupCounter(v, strand, max_del)

      # skip ambiguous or soft-masked sequences
      if current_rec.base not in ref_nts:
          continue

      if n_good_variants == 0:
          n_good_variants += 1
          ic = IndelCache(current_rec.chrom)
          previous_rec = current_rec 
          continue

      if ic.chrom != current_rec.chrom:
          if debug:
            print("after processing {},  {} records remain".format(ic.chrom, 
                        ic.n()))
            print(ic)
          ic.clear()
          ic = IndelCache(current_rec.chrom)
      
      if current_rec.has_indel:
          ic.add(current_rec)

      if current_rec.pos in ic.d:
         icounts = ic.get_indel_depth(current_rec)
         current_rec.indelcount += icounts
     
      if current_rec.pos in ic.refdepth:
         ic_counts = ic.get_indel_skipped_coverage(current_rec)
         current_rec.indelcovcount += ic_counts
      
      all_bases = previous_rec
      previous_rec = current_rec 
      
      if all_bases.get_depth() >= min_depth:
          out_fh.write(str(all_bases))
  
  if previous_rec is not None:
      if previous_rec.get_depth() >= min_depth:
          out_fh.write(str(previous_rec))
  
  if debug:
    print("after processing {},  {} and {} records remain".format(ic.chrom, 
                ic.n()[0], ic.n()[1]))
    print(ic)
  
  out_fh.close()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="""
    description""")

    parser.add_argument('-p',
                        '--pileup',
                        help ="""pileup file input, use "-" for stdin
                        \n""",
                        required = True)
                        
    parser.add_argument('-o',
                        '--output',
                        help="""output file name
                        \n""",
                        required = True,
                        default = "counts.tsv.gz")
    
    parser.add_argument('-d',
                        '--depth',
                        help="""minimum depth to report
                        \n""",
                        required = False,
                        default = 1,
                        type = int)
    
    parser.add_argument('-r',
                        '--complement',
                        help = """
                        if -r is set then the reported nucleotides will be complemented. 
                        (default: %(default)s)
                        \n""",
                        action = 'store_true')

    parser.add_argument('-l',
                        '--deletion_length',
                        help = """
                        do not count deletions greater than or equal
                        to deletion_length (default: %(default)s)
                        \n""",
                        required = False,
                        default = 4, 
                        type = int)

    parser.add_argument('--verbose',
                        help = """
                        (default: %(default)s)
                        \n""",
                        action = 'store_true')
    
    args = parser.parse_args()

    if args.pileup == "-":
        fh = sys.stdin
    else:
        fh = open(args.pileup, 'r')

    if args.output == "-":
        fo = sys.stdout
    else:
        fo = gzip.open(args.output, 'wt')

    mpileup_to_counts(fh, fo,
            int(args.depth), args.complement,
            debug = args.verbose, max_del = args.deletion_length,
            header = False)

