#!/usr/bin/env python3
import argparse
import os 
import sys
import pysam
import pdb
import re

from cyvcf2 import VCF
from collections import Counter
from yaml import load, dump

from utils import * 

import pdb

bp_dict = {
        "A" : "T",
        "T" : "A",
        "C" : "G",
        "G" : "C"}

def variant_to_pileup(ref, variants):
    var_dict = {}
    for var in variants:
        if var == '<*>':
            continue
        len_diff = len(ref) - len(var)
        if len_diff == 0:
            sign = ""
        elif len_diff > 0:
            sign = "-"
        else:
            sign = "+"

        res = ref[0] + sign
        if(len_diff > 0):
          res += str(len_diff) + ('N' * len_diff)
        elif len_diff < 0:
            res += str(abs(len_diff)) + var[len_diff:]
        var_dict[var] = res
    return(var_dict)

def pileup_variant_to_vcf(ref, pileup_variants, vcf_variants):
    """ just handle indels for now""" 
    var_dict = {}
    offset_regex = re.compile('[0-9]+')
    for var in pileup_variants:
        if var.startswith("^"):
            continue

        if "-" in var:
            try:
                offset = int(offset_regex.search(var).group())
            except AttributeError:
                sys.exit("unable to parse variant " + var)
            res = ref[0] + ref[1 + offset:]
        elif "+" in var:
            try:
                offset = offset_regex.search(var).end()
            except AttributeError:
                sys.exit("unable to parse variant " + var)
            res = ref[0] + var[offset:] + ref[1:]
        else:
            continue

        if res not in vcf_variants:
            res = var
       
        var_dict[var] = res
    return(var_dict)

        
def count_variants(pileup, variant, warn_allele_mismatch = True):
    """
    given pysam alignments at a position
    determine how many have the variant
    variant is supplied as cyvcf2 variant class
    (useful for counting indels with are not counted correctly by
    bcftools)
    """
    if (len(variant.ALT) == 1) and (variant.ALT[0] == '<*>'):
        return
        
    n = 0
    for pileupcolumn in pileup :
        if n > 1:
            sys.exit("only a single pileup position should be queried")
        n += 1

        pup = pileupcolumn.get_query_sequences(mark_matches =
                              True,mark_ends = True,
                              add_indels = True)

        pup = [x.upper() for x in pup]
        
        frequency = Counter()
        for val in pup:
            if val.startswith("^"):
                continue
            if "-" in val: 
                frequency[val] += 1
            elif "+" in val:
                # swap first nucleotide with ref, samtools will report as
                # snp if not the same
                val = val.replace(val[val.index("+") - 1], variant.REF[0], 1)
                frequency[val] += 1
        
        vcf_variant_dict = pileup_variant_to_vcf(variant.REF,
                                                 list(frequency.keys()),
                                                 variant.ALT)
        res = Counter()
        for k,v in frequency.items():
            res[vcf_variant_dict[k]] = v
        
        if warn_allele_mismatch:
          alts = set(variant.ALT)
          new_vars = []
          for x in res.keys():
            if x not in alts:
                new_vars.append(x)
            else:
                alts.remove(x)
          if len(alts) > 0:
              print(variant.CHROM, 
                    variant.POS,
                    "not all variant were found: ", 
                    ",".join(alts),
                    file = sys.stderr)
          if len(new_vars) > 0:
              print(variant.CHROM,
                    variant.POS,
                    "the following new variants were found: ",
                    ",".join(new_vars),
                    file = sys.stderr)
        
    return(res)


def scan_bam(bam_fn, vcf_fn, only_indels, pup_args): 
    """ returns fileobject to pileup output 
    args:
        bam = path to bam
    """
    
    bam_fh = pysam.AlignmentFile(bam_fn, "rb")
    for variant in VCF(vcf_fn):

      if only_indels and variant.INFO.get('INDEL') is None:
          continue
    
      pileup_itr = bam_fh.pileup(variant.CHROM, variant.start,
              variant.start + 1, 
              truncate = True, 
              stepper = 'all',
              **pup_args)
      res = count_variants(pileup_itr, variant)    
      
      for k,v in res.items():
          print(variant.REF, k, v)
    bam_fh.close()

def main():
    
    parser = argparse.ArgumentParser(description="""
    count reads with specific variants
    """,
    formatter_class = argparse.RawTextHelpFormatter )

    parser.add_argument('-b',
                        '--bam',
                        help ="""indexed bam file input,
                        by default the bam will be split into
                        forward and reverse bams,
                        pass the forward and reverse bams
                        to bypass splitting the bams
                        (i.e. forward.bam,reverse.bam)
                        \n""",
                        required = True)
    parser.add_argument('-v',
                        '--vcf',
                        help = """\
                        VCF file with variants of interest
                        \n""",
                        required = True)
    parser.add_argument('-i',
                        '--indels',
                        help = """\
                        only process indels
                        \n""",
                        action = 'store_true')
    parser.add_argument('-p',
                        '--pileup-args',
                        help = """\
                        a config file (pileup.yaml) in
                        the src/ directory is used to set 
                        default pileup arguments. 
                        Provide path to custom yaml file here 
                        or edit pileup.yaml to change defaults. 
                        \n""",
                        required = False)
 
        
    args = parser.parse_args()

    if args.pileup_args:
      config_fn = args.pileup_args
    else:
      config_fn = os.path.join(os.path.dirname(os.path.abspath(__file__)),
              "pileup.yaml")
    
    with open(config_fn) as f:
        config_data = f.read()

    pup_args = load(config_data)
    
    bam_name = scan_bam(args.bam, args.vcf, args.indels, pup_args)
                
if __name__ == '__main__': main()

