#!/usr/bin/env python3
import argparse
import textwrap
import os 
import sys
import pysam
import pandas as pd
import numpy as np

import norm_methods
import utils

def parse_bed(fh):

    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.split()
        
        if len(fields) < 3:
            continue 

        fields[1] = int(fields[1])
        fields[2] = int(fields[2])
        yield fields

bpmap = {
        "A" : "T",
        "T" : "A",
        "G" : "C",
        "C" : "G",
        "N" : "N"}

class Mapfile:
    
    def __init__(self, nt):
        self.mm = -999.0
        self.se = 0.0
        self.nt = nt
    
    def __str__(self):
        return "{}\t{}\t{}".format(self.mm, self.se, self.nt)
    
    def to_list(self):
        return [self.mm, self.se, self.nt]

class Structinfo:

    def __init__(self, string, m_col, se_col):
        vals = string.split()
        self.chr = vals[0]
        self.pos = int(vals[1])
        self.strand = vals[2]
        self.ref_base = vals[3]
        self.depth = int(vals[4])
        self.refcount = int(vals[5])
        self.mm = float(vals[m_col])
        self.se = float(vals[se_col])
    
    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}\t{}".format(self.chr, 
                                        self.pos,
                                        self.ref_base,
                                        self.strand,
                                        self.mm,
                                        self.se)

class Structpileup:

    def __init__(self, string):
        vals = string.split()
        self.chr = vals[0]
        self.pos = int(vals[1])
        self.strand = vals[2]
        self.ref_base = vals[3]
        self.depth = int(vals[4])
        self.refcount = int(vals[5])
        self.acount = int(vals[6])
        self.ccount = int(vals[7])
        self.gcount = int(vals[8])
        self.tcount = int(vals[9])
        self.ncount = int(vals[10])
        self.mmcount = int(vals[11])
        self.delcount = int(vals[12])
        self.inscount = int(vals[13])
        self.mmratio = float(vals[14])
        self.insratio = float(vals[15])
        self.delratio = float(vals[16])
        self.mm = float(vals[23])
        self.se = float(vals[24])
 
    def __str__(self):
        return "\t".join([str(x) for x in vars(self).values()])

def get_global_norm(gene_fn, tbx_fo, m_col, se_col):
    mutation_rates = []
    for tx in parse_bed(gene_fn):
        
        has_strand = False
        if len(tx) >= 6:
            if tx[5] in ["+", "-"]:
                has_strand = True
            else:
                print("unrecognized character in column 6, expecting + or - \n ignoring strand", file = sys.stderr)
        
        for row in tbx_fo.fetch(tx[0], tx[1], tx[2]):
            
            vals = Structinfo(row, m_col, se_col)

            if has_strand and vals.strand != tx[5]:
              continue
          
            # ignore - pileup features if no strandedness requested from
            # interval
            if not has_strand and vals.strand == "-":
              continue

            mutation_rates.append(vals.mm)
    mutation_rates = np.array(mutation_rates)
    norm_factor = norm_methods.calc_global_norm(mutation_rates)
    
    gene_fn.seek(0)
    return norm_factor

def main():
    
    parser = argparse.ArgumentParser(description="""
        Convert genomic coverage values from tabix table to .map file for
        SuperFold, constraints file for RNAfold, or 
        """,
        formatter_class = argparse.RawTextHelpFormatter )

    parser.add_argument('-t',
                        '--tabix-file',
                        help ="""tabix indexed pileup file
                        \n""",
                        required = True)

    parser.add_argument('-g',
                        '--gene_bed',
                        help = textwrap.dedent("""\
                        bedfile with transcript coordinates to convert to
                        map (zero-based intervals)
                        \n"""),
                        required = True)
    
    parser.add_argument('-f',
                        '--fasta',
                        help = textwrap.dedent("""\
                        fasta file, needed to fill in nucleotides from
                        regions without enough depth to be present in
                        pileup table
                        \n"""),
                        required = True)
    parser.add_argument('-m',
                        '--mutation_col',
                        help = textwrap.dedent("""\
                        column index of shape reactivities (1-based)
                        \n"""),
                        required = False,
                        default = 24, 
                        type = int)
    parser.add_argument('-e',
                        '--se_col',
                        help = textwrap.dedent("""\
                        column index of shape reactivity standard error (1-based)
                        \n"""),
                        required = False,
                        default = 25,
                        type = int)
    
    parser.add_argument('-T',
                        '--outtype',
                        help = textwrap.dedent("""\
                        Output type, either (M) map,(S) shape,or (C) constraint
                        \n"""),
                        required = False,
                        default = "M")

    parser.add_argument('-n',
                        '--norm',
                        help = textwrap.dedent("""\
                        Normalized data. Windsorize for constraint,
                        or boxplot norm for shape and map.
                        \n"""),
                        action = "store_true")

    parser.add_argument('-gn',
                        '--global_norm',
                        help = textwrap.dedent("""\
                        Estimate normalization cutoff using all supplied
                        bed interval regions (for shape and Map only)
                        \n"""),
                        action = "store_true")
    parser.add_argument('-nc',
                        '--norm_cutoff',
                        help = textwrap.dedent("""\
                        Cutoff value to determine paired or unpaired for
                        constraint or normalization threshold for MAP or
                        shape values
                        \n"""),
                        required = False,
                        type = float)

    parser.add_argument('-o',
                        '--outdir',
                        help = textwrap.dedent("""\
                        output directory for files
                        \n"""),
                        required = False,
                        default = "./")

    args = parser.parse_args() 
    
    gene_fn = open(args.gene_bed)
    
    tbx = pysam.TabixFile(args.tabix_file)
    output_dir = args.outdir
    
    m_col = int(args.mutation_col) - 1
    se_col = int(args.se_col) - 1

    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
   
    fa_fh = pysam.FastaFile(args.fasta)
    
    norm_data = args.norm
    norm_cutoff = args.norm_cutoff
    outtype = args.outtype
        
    if outtype in ["M", "S"]:
        if norm_data and not norm_cutoff:
            norm_cutoff = get_global_norm(gene_fn, tbx, m_col, se_col)
     
    for tx in parse_bed(gene_fn):
        
        d = {}

        #get nucleotides from fasta
        seqs = fa_fh.fetch(tx[0], tx[1], tx[2])
        seqs = seqs.upper()
        
        has_strand = False
        if len(tx) >= 6:
            if tx[5] in ["+", "-"]:
                has_strand = True
            else:
                print("unrecognized character in column 6, expecting + or - \n ignoring strand", file = sys.stderr)
        
        revcomp = False
        if has_strand and tx[5] == "-":
            seqs = [bpmap[x] for x in seqs]
            revcomp = True

        for idx,nt in enumerate(seqs, 1):
            d[idx] = Mapfile(nt)
        
        for row in tbx.fetch(tx[0], tx[1], tx[2]):
            
            vals = Structinfo(row, m_col, se_col)

            if has_strand and vals.strand != tx[5]:
              continue
          
            # ignore - pileup features if no strandedness requested from
            # interval
            if not has_strand and vals.strand == "-":
              continue

            new_start = max(vals.pos - tx[1], 1) 

            d[new_start].mm = vals.mm
            d[new_start].se = vals.se
            d[new_start].nt = vals.ref_base
        
        res_list = []
        if revcomp:
            vals = list(d.values())[::-1]
            for idx,val in enumerate(vals):
                nidx = idx + 1
                res_list.append([nidx] +  val.to_list())
        else:
            for key,val in d.items():
                res_list.append([key] + val.to_list())

        df = pd.DataFrame(res_list)
        

        if norm_data and outtype == "C":
            nvals = df.shape[0]
            nas = np.where(df[1] == -999)[0]
            vals = np.where(df[1] != -999)[0]

            res = norm_methods.windsorize(df[1][vals]).tolist()
             
            norm_react = []
            for i in range(nvals):
              if i in nas:
                  norm_react.append(-999.0)
              elif i in vals:
                  new_val = res.pop(0)
                  norm_react.append(new_val)
              else:
                  sys.exit("shouldn't have reached this point")
            df["norm_vals"] = norm_react 
            
            if norm_cutoff:
              df["constraints"] = norm_methods.get_constraints(df["norm_vals"], simple_cutoff = norm_cutoff)             
            else:
              df["constraints"] = norm_methods.get_constraints(df["norm_vals"]) 

        if outtype == "M":
            fout = open(os.path.join(output_dir, tx[3] + "_coverage.map"), 'w')
            for line in res_list:
                if norm_data:
                    if line[1] == -999:
                        pass
                    else:
                        line[1] = float(line[1]) / norm_cutoff
                        line[2] = float(line[2]) / norm_cutoff
                line = [str(x) for x in line]
                fout.write("\t".join(line) + "\n")
            fout.close()
        elif outtype == "S" :
            pass

        elif outtype == "C":
          utils.to_fasta(df, 
                os.path.join(output_dir, tx[3] + ".constraints"),
                nuc_col_idx = 3,
                header_str = tx[3],
                other_data_idx = 5)

if __name__ == '__main__': main()


