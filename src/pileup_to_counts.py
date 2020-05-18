#! /usr/bin/env python3
import pysam
import numpy as np
import argparse
import sys
import os
import sys
import gzip
from cyvcf2 import VCF
from collections import Counter

import variant_counter
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
  
    alleles_to_count = ["A","T","C","G","indel", "indel_cov"]
    alleles_dict = Counter()
    
    for base in alleles_to_count:
        alleles_dict[base] = 0  
    
    return alleles_dict

def is_snp(v):
    "specific for bcftools mpileup vcfs"
    if len(v.REF) > 1: return False
    for i in range(0, len(v.ALT)):
        if not v.ALT[i] in ("A", "C", "G", "T", "<*>"):
            return False
        return (len(v.REF) + len(v.ALT)) > 1

def parse_pileup_variant(v):
    """
    return deletion length (or 1 for insertions)
    from novel variants discovered during pileup
    """
    
    offset_regex = re.compile('[0-9]+')
    if "+" in v:
        offset = -1 
    elif "-" in v:
        offset = int(offset_regex.search(v).group())
    else:
        sys.exit(v, " is not a deletion or insertion")
    return offset

class IndelCache():
    """cache to hold vcf records for indels. Used to ensure
    that indel position reporting is left aligned if forward stranded
    alignments were used and right aligned if reverse stranded.
    Thinking that RT will hit an adduct from 3'-> 5'. For ambigious
    indels, there isn't a way to precisley localize, however
    left alignment is used in some publications"""

    def __init__(self, chrom = None):
        self.d = Counter()
        self.refdepth = Counter()
        self.chrom = chrom
        
    def add(self, ac):
        
        if ac.multiallelic:
            for var, depth in ac.new_alleles.items():
              
              if var not in ac.v.ALT:
                  del_length = parse_pileup_variant(var)
              else:    
                  del_length = len(ac.v.REF) - len(var)

              if ac.strand == "+":
                  if del_length < 0:
                      # insertion
                      # annotate to position 5' of insertion site
                      self.d[ac.v.POS] += depth
                  else:
                      # deletion
                      # annotation to 3' most position of deletion
                      self.d[ac.v.POS + del_length] += depth
                      # add to coverage
                      for i in range(ac.v.POS + 1, ac.v.POS + 1 + del_length):
                          self.refdepth[i] += depth
              else:
                  if del_length < 0:
                      # insertion
                      # annotate to position 3' of insertion 
                      self.d[ac.v.POS + 1] += depth
                  else:
                      # deletion
                      # annotate to 5' most position of deletion (5' in
                      # vcf)
                      self.d[ac.v.POS + 1] += depth
                      for i in range(ac.v.POS + 1, ac.v.POS + 1 + del_length): 
                          self.refdepth[i] += depth

        else:
          
          del_length = len(ac.v.REF) - len(ac.v.ALT[0])
          
          if ac.strand == "+":
              # report left aligned position, offset by length of deletion
              if del_length < 0:
                  # insertion
                  self.d[ac.v.POS] += ac.allele_counts["indel"] 
              else:
                  # deletion
                  self.d[ac.v.POS + del_length] += ac.allele_counts["indel"]
                  for i in range(ac.v.POS + 1, ac.v.POS + 1 + del_length):
                      self.refdepth[i] += ac.allele_counts["indel"] 
          else:
              if del_length < 0:
                  # insertion
                  # annotate to position 3' of insertion 
                  self.d[ac.v.POS + 1] += ac.allele_counts["indel"]
              else:
                  # deletion
                  # annotate to 5' most position of deletion (5' in
                  # vcf)
                  self.d[ac.v.POS + 1] += ac.allele_counts["indel"]
                  for i in range(ac.v.POS + 1, ac.v.POS + 1 + del_length): 
                      self.refdepth[i] += ac.allele_counts["indel"]
              

    def n(self):
        return len(self.d), len(self.refdepth)

    def get_indel_depth(self, ac):
        return self.d.pop(ac.v.POS, 0)
    
    def get_indel_skipped_coverage(self, ac):
        return self.refdepth.pop(ac.v.POS, 0)
    
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

class AlleleCounter():
    "Single sample allele counter for bcftools mpileup VCF"
    
    def __init__(self, variant, bam, pup_args, strand = "+"):
        self.v = variant
        self.b = bam
        self.pileup_args = pup_args
        self.strand = strand
        self.ref_depths = np.array(variant.gt_ref_depths[0])
        self.alt_depths = np.array(variant.format('AD')[0][1:])
        self.allele_counts = allele_counter()
        self.is_indel = False
        self.indel_pos = None
        self.multiallelic = False
        self.new_alleles = None
        self.get_allele_counts()
        
    def __str__(self):
        return "\t".join([str(x) for x in self.to_counts()]) + "\n"
    
    def depth(self):
        # return counts excluding indel values but including indel coverage
        # values
        return sum(self.allele_counts.values()) - self.allele_counts["indel"]

    def to_counts(self):
        depth = self.depth() 
        indel_counts = self.allele_counts["indel"]
        ref_counts = self.allele_counts[self.v.REF]
        mm_counts = self.get_mismatch_counts(self.v.REF)

        if self.strand == "-":
            ref_nt = bp_comp[self.v.REF]
            
            out = [self.v.CHROM,
            self.v.POS,
            self.strand,
            ref_nt,
            depth,
            ref_counts,
            self.allele_counts[bp_comp["A"]],
            self.allele_counts[bp_comp["C"]],
            self.allele_counts[bp_comp["G"]],
            self.allele_counts[bp_comp["T"]],
            mm_counts,
            indel_counts
            ]
        else:    
            out = [self.v.CHROM,
            self.v.POS,
            self.strand,
            self.v.REF,
            depth,
            ref_counts,
            self.allele_counts["A"],
            self.allele_counts["C"],
            self.allele_counts["G"],
            self.allele_counts["T"],
            mm_counts,
            indel_counts
            ]
        return out
    
    def query_indels(self):
        # need to query pileup using pysam
        pileup_itr = self.b.pileup(self.v.CHROM, 
                                   self.v.start,
                                   self.v.start + 1,
                                   truncate = True,
                                   **self.pileup_args)
        
        return variant_counter.count_variants(pileup_itr, self.v)
        

    def get_allele_counts(self):
        if self.v.INFO.get("INDEL") is not None:
            # the AD field is not correct for indels :(
            # See open issue
            # https://github.com/samtools/bcftools/issues/912
            # IDV is correct for single indels, but need to reparse pileup
            # for multiple indels
            if len(self.v.ALT) > 1:
                indel_depths = self.query_indels()
                self.multiallelic = True
                self.new_alleles = indel_depths
                self.allele_counts["indel"] += sum(indel_depths.values())
            else:    
                self.allele_counts["indel"] += self.v.INFO.get("IDV")
            self.is_indel = True
            # set indel position based on strandedness
            # if + just set to next position 
            # if - set to length of REF 
            # need to look into if right alignment of the indel is
            # going to mess this up

            self.indel_pos = self.v.POS + 1
            
        elif is_snp(self.v):
            self.allele_counts[self.v.REF] += self.ref_depths
            for idx, snp in enumerate(self.v.ALT):
                if snp == "<*>":
                    continue
                self.allele_counts[self.v.ALT[idx]] += self.alt_depths[idx]
    
    def get_mismatch_counts(self, ref_nt):
        nt = {"A", "T", "C", "G"}
        if ref_nt in nt:
            nt.remove(ref_nt)
        return sum([self.allele_counts[x] for x in nt])
        
    def add_indel_depth(self, new_variant):
        new_counter = AlleleCounter(new_variant, self.strand)
        new_counter.allele_counts["indel"] = self.allele_counts["indel"]
        new_counter.allele_counts["indel_cov"] = self.allele_counts["indel_cov"]

        new_counter.ref_depths += self.ref_depths
        return new_counter

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
     "indelcount"]
  
  fo.write("\t".join(tbl_cols) + "\n")

def vcf_to_counts(vcf_fn, bam_fn, out_fn, pileup_args, min_depth = 10, return_comp = False,
        n_records= -1, debug = True, region = None): 
      
  if return_comp:
      strand = "-"
  else:
      strand = "+"
  
  if out_fn == "-":
    fo = sys.stdout
  else:
    fo = gzip.open(out_fn, 'wt')
  
  write_header(fo)
  vcf = VCF(vcf_fn)
  bam = pysam.AlignmentFile(bam_fn)
  n_good_variants = 0

  ic = IndelCache()
  previous_rec = None
  for i,v in enumerate(vcf(region)):

      if i == n_records:
          break

      if n_good_variants == 0:
          n_good_variants += 1
          ic = IndelCache(v.CHROM)
          previous_rec = AlleleCounter(v, bam, pileup_args, strand)
          continue
      
      if ic.chrom != v.CHROM:
          if debug:
            print("after processing {},  {} records remain".format(ic.chrom, 
                        ic.n()))
            print(ic)
          ic.clear()
          ic = IndelCache(v.CHROM)
      
      if v.INFO.get('DP') < min_depth and v.POS not in ic.d and v.POS not in ic.refdepth:
          continue
      #pdb.set_trace()
      current_rec = AlleleCounter(v, bam, pileup_args, strand)
      
      if current_rec.is_indel:
          ic.add(current_rec)

          if previous_rec.v.POS in ic.d:
              # only happens for insertions
              ins_count = ic.d.pop(previous_rec.v.POS)
              previous_rec.allele_counts["indel"] += ins_count
          all_bases = previous_rec
          # dont report indel as single record
          previous_rec = current_rec
      else:
          if current_rec.v.POS in ic.d:
              # add indel depths to current record
              # pop POS from IndelCache()
              icounts = ic.get_indel_depth(current_rec)
              current_rec.allele_counts["indel"] += icounts
              if current_rec.v.POS in ic.refdepth:
                ic_counts = ic.get_indel_skipped_coverage(current_rec)
                current_rec.allele_counts["indel_cov"] += ic_counts
              all_bases = previous_rec
              previous_rec = current_rec
          else:
              if current_rec.v.POS in ic.refdepth:
                  icounts = ic.get_indel_skipped_coverage(current_rec)
                  current_rec.allele_counts["indel_cov"] += icounts
              if previous_rec.is_indel:
                  previous_rec = current_rec
                  continue
              all_bases = previous_rec
              previous_rec = current_rec 
      
      if all_bases.is_indel:
            continue
      
      if all_bases.depth() >= min_depth:
          fo.write(str(all_bases))
  if previous_rec is not None and not previous_rec.is_indel:
      if previous_rec.depth() >= min_depth:
          fo.write(str(previous_rec))
  
  if debug:
    print("after processing {},  {} and {} records remain".format(ic.chrom, 
                ic.n()[0], ic.n()[1]))
    print(ic)
  
  vcf.close()
  fo.close()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="""
    description""")

    parser.add_argument('-v',
                        '--vcf',
                        help ="""vcf file input, use "-" for stdin
                        \n""",
                        required = True)
    
    parser.add_argument('-b',
                        '--bam',
                        help ="""indexed bam file input
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
    parser.add_argument('-q',
                        '--query-region',
                        help = """
                        query supplied region e.g. ch17:38544980-38544982
                        requires a vcf indexed by bcftools
                        """,
                        required = False)
    
    parser.add_argument("--pileup-args",
                      nargs='+',
                      action=utils.kvdictAppendAction,
                      metavar="KEY=VALUE",
                      help="""Add key/value pileup arguments to overwrite
                      defaults, also can use --pileup-arg-file to specify
                      args""")
    parser.add_argument("--pileup-arg-fn",
                      help="""
                      specify custom pileup.yaml file to overwrite default
                      arguments
                      """,
                      required = False)
    parser.add_argument('--verbose',
                        help = """
                        (default: %(default)s)
                        \n""",
                        action = 'store_true')
    
    args = parser.parse_args()
    pileup_arguments = utils.get_pileup_args(args.pileup_arg_fn, args.pileup_args)
    nrecords = -1 

    vcf_to_counts(args.vcf, args.bam, args.output,
            pileup_arguments, int(args.depth), args.complement,
            debug = args.verbose, region = args.query_region, n_records =
            nrecords)

