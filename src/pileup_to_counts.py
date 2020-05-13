#! /usr/bin/env python3
import pysam
from cyvcf2 import VCF
import numpy as np
import argparse
import sys
import os
import gzip
import sys
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
        offset = 1 
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

    def __init__(self, bam, chrom = None):
        self.d = Counter()
        self.refdepth = Counter()
        self.chrom = chrom
        self.b = bam
        
    def add(self, ac):
        # affected nucleotide is first deleted or first insetion pos from
        # 3' end
        if ac.multiallelic:
            for var, depth in ac.new_alleles.items():
              if ac.strand == "+":
                  if var not in ac.v.ALT:
                      del_length = parse_pileup_variant(var)
                  else:    
                  # report left aligned position, offset by length of deletion
                  # or report insertion at position immediately 5'
                      del_length = max(len(ac.v.REF) - len(var), 1)
                  self.d[ac.v.POS + del_length] += depth
                  for i in range(ac.v.POS + 1, ac.v.POS + del_length):
                      self.refdepth[i] += depth
              else:
                  if var not in ac.v.ALT:
                      del_length = parse_pileup_variant(var)
                  # report right aligned position
                  else:
                      del_length = max(len(var) - 1, 1)
                  right_align_pos = ac.v.POS + del_length 
                  self.d[right_align_pos] += depth

                  for i in range(ac.v.POS + 1, right_align_pos): 
                      self.refdepth[i] += depth

        else:
          if ac.strand == "+":
              # report left aligned position, offset by length of deletion
              del_length = max(len(ac.v.REF) - len(ac.v.ALT[0]), 1)
              self.d[ac.v.POS + del_length] += ac.allele_counts["indel"] 
              for i in range(ac.v.POS + 1, ac.v.POS + del_length):
                 self.refdepth[i] += ac.allele_counts["indel"] 
          else:
              # report right aligned position
              right_align_pos = ac.v.POS + max(len(ac.v.ALT[0]) - 1, 1)
              self.d[right_align_pos] += ac.allele_counts["indel"] 
              
              for i in range(ac.v.POS + 1, right_align_pos): 
                  self.refdepth[i] += ac.allele_counts["indel"]

    def n(self):
        return len(self.d), len(self.refdepth)

    def get_indel_depth(self, ac):
        _ = self.get_indel_skipped_coverage(ac)
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
    
    def __init__(self, variant, bam, strand = "+"):
        self.v = variant
        self.b = bam
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
    
    def to_counts(self):
        depth = sum(self.allele_counts.values()) 
        indel_counts = self.allele_counts["indel"]
        ref_counts = self.allele_counts[self.v.REF]
        mm_counts = depth - indel_counts - ref_counts - self.allele_counts["indel_cov"]
        
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
                                   **utils.pileup_args)
        
        return variant_counter.count_variants(pileup_itr, self.v)
        

    def get_allele_counts(self):
        if self.v.INFO.get("INDEL") is not None:
            # the AD field is not correct for indels :(
            # See open issue
            # https://github.com/samtools/bcftools/issues/912
            
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

def vcf_to_counts(vcf_fn, bam_fn, out_fn, min_depth = 10, return_comp = False,
        n_records= -1, debug = True, region = None): 
      
  if return_comp:
      strand = "-"
  else:
      strand = "+"
  
  fo = gzip.open(out_fn, 'wt')
  
  write_header(fo)
  vcf = VCF(vcf_fn)
  bam = pysam.AlignmentFile(bam_fn)
  n_good_variants = 0

  ic = IndelCache(bam)
  previous_rec = None
  for i,v in enumerate(vcf(region)):

      if i == n_records:
          break
      
      if n_good_variants == 0:
          n_good_variants += 1
          ic = IndelCache(bam, v.CHROM)
          previous_rec = AlleleCounter(v, bam, strand)
          continue
      
      if ic.chrom != v.CHROM:
          if debug:
            print("after processing {},  {} records remain".format(ic.chrom, 
                        ic.n()))
            print(ic)
          ic.clear()
          ic = IndelCache(bam, v.CHROM)

      current_rec = AlleleCounter(v, bam, strand)
      
      if current_rec.is_indel:
        #  pdb.set_trace()
          ic.add(current_rec)
          all_bases = previous_rec
          # dont report indel as single record
          previous_rec = current_rec
      else:
          if current_rec.v.POS in ic.d:
              # add indel depths to current record
              # pop POS from IndelCache()
              icounts = ic.get_indel_depth(current_rec)
              current_rec.allele_counts["indel"] += icounts
              # need to consume reference depth due to right alignment of
              # indel, unless next to indel
#              if not previous_rec.is_indel:
#                current_rec.allele_counts[current_rec.v.REF] -= icounts
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
      # need to check if indel cache is still present at end of contig
      # debug with printing

      if sum(all_bases.allele_counts.values()) >= min_depth:
          fo.write(str(all_bases))
  
  if previous_rec is not None and not previous_rec.is_indel:
      if sum(all_bases.allele_counts.values()) >= min_depth:
          fo.write(str(previous_rec))
  
  if debug:
    print("after processing {},  {} and {} records remain".format(ic.chrom, 
                ic.n()[0], ic.n()[1]))
    print(ic)
  
  vcf.close()
  fo.close()

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="""
    description""",
    formatter_class = argparse.RawTextHelpFormatter )

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
    parser.add_argument('--verbose',
                        help = """
                        (default: %(default)s)
                        \n""",
                        action = 'store_true')
    
    args = parser.parse_args()
    
    nrecords = -1 
    vcf_to_counts(args.vcf, args.bam, args.output, int(args.depth), args.complement,
            debug = args.verbose, region = args.query_region, n_records = nrecords)

