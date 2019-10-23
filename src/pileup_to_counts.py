#! /usr/bin/env python3


from cyvcf2 import VCF
import numpy as np
import argparse
import sys
import os
import gzip
import sys
from collections import Counter

bp_comp = {
        "A" : "T",
        "T" : "A",
        "G" : "C",
        "C" : "G",
        "N" : "N"
        }

ref_nts = ["A", "T", "C", "G"]

def allele_counter():
  
    alleles_to_count = ["A","T","C","G","indel"]
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
    
class AlleleCounter():
    "Single sample allele counter for bcftools mpileup VCF"
    
    def __init__(self, variant, strand = "+"):
        self.v = variant
        self.strand = strand
        self.ref_depths = np.array(variant.gt_ref_depths[0])
        self.alt_depths = np.array(variant.format('AD')[0][1:])
        self.allele_counts = allele_counter()
        self.is_indel = False
        self.indel_pos = None
        self.get_allele_counts()
        
    def __str__(self):
        return "\t".join([str(x) for x in self.to_counts()]) + "\n"
    
    def to_counts(self):
        depth = sum(self.allele_counts.values())
        indel_counts = self.allele_counts["indel"]
        ref_counts = self.allele_counts[self.v.REF]
        mm_counts = depth - indel_counts - ref_counts
        
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
    
    def get_allele_counts(self):
        if self.v.INFO.get("INDEL") is not None:
            # the AD field is not correct for indels :(
            # See open issue
            # https://github.com/samtools/bcftools/issues/912
            self.allele_counts["indel"] += self.v.INFO.get("IDV")
            self.is_indel = True
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

def vcf_to_counts(vcf_fn, out_fn, min_depth = 10, return_comp = False, n_records= -1): 
      
  if return_comp:
      strand = "-"
  else:
      strand = "+"
  
  fo = gzip.open(out_fn, 'wt')
  
  write_header(fo)
  vcf = VCF(vcf_fn)
  n_good_variants = 0
  for i,v in enumerate(vcf):

      if i == n_records:
          break
      
      if n_good_variants == 0:
          n_good_variants += 1
          ac = AlleleCounter(v, strand)
          continue
      
      if (v.CHROM, v.POS) == (ac.v.CHROM, ac.v.POS):
          # check whether previous nt was a an indel
          # if so, add indel counts 
          tmp_ac = AlleleCounter(v, strand)
          if tmp_ac.is_indel:
              all_bases = ac
              ac = tmp_ac
          else:
              sys.exit("not sure what to do")
      else:
          if ac.is_indel and (v.CHROM == ac.v.CHROM and v.POS == ac.indel_pos):
              ac = ac.add_indel_depth(v)
              continue
          else:
              all_bases = ac
              ac = AlleleCounter(v, strand)  
       
      if sum(all_bases.allele_counts.values()) >= min_depth:
          fo.write(str(all_bases))
      
  fo.write(str(ac))
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
                        
    parser.add_argument('-o',
                        '--output',
                        help="""output file name
                        \n""",
                        required = True,
                        default = "counts.tsv.gz")
    
    parser.add_argument('-d',
                        '--depth',
                        help="""minimum depth tp report
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
    
    args = parser.parse_args()
    
    
    vcf_to_counts(args.vcf, args.output, int(args.depth), args.complement)

