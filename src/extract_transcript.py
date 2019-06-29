#!/usr/bin/env python3
import argparse
import textwrap
import os 
import sys
import pysam

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

class Structinfo:

    def __init__(self, string):
        vals = string.split()
        self.chr = vals[0]
        self.pos = int(vals[1])
        self.strand = vals[2]
        self.ref_base = vals[3]
        self.depth = int(vals[4])
        self.refcount = int(vals[5])
        self.mm = float(vals[17])
        self.se = float(vals[18])
    
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

def main():
    
    parser = argparse.ArgumentParser(description="""
        Convert genomic coverage values from genomic bedgraph to .map file for Shape
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

    parser.add_argument('-o',
                        '--outdir',
                        help = textwrap.dedent("""\
                        output directory for map files
                        \n"""),
                        required = False,
                        default = "./")

    args = parser.parse_args() 
    
    gene_fn = open(args.gene_bed)
    
    tbx = pysam.TabixFile(args.tabix_file)
    output_dir = args.outdir

    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
   
    fa_fh = pysam.FastaFile(args.fasta)
    
    for tx in parse_bed(gene_fn):
        fout = open(os.path.join(output_dir, tx[3] + "_coverage.map"), 'w')
        
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
            
            vals = Structinfo(row)

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

        if revcomp:
            vals = list(d.values())[::-1]
            for idx,val in enumerate(vals):
                nidx = idx + 1
                fout.write(str(nidx) + "\t" + str(val) + "\n")
        else:
            for key,val in d.items():
                fout.write(str(key) + "\t" + str(val) + "\n")
 
if __name__ == '__main__': main()

