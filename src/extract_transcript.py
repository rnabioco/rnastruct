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

        fields = fields[0:4]
        fields[1] = int(fields[1])
        fields[2] = int(fields[2])
        yield fields

def main():
    
    parser = argparse.ArgumentParser(description="""
    Convert genomic coverage values from genomic bedgraph to .map file for Shape
    """,
    formatter_class = argparse.RawTextHelpFormatter )

    parser.add_argument('-b',
                        '--bedgraph',
                        help ="""bedgraph
                        \n""",
                        required = True)
    parser.add_argument('-g',
                        '--gene_bed',
                        help = textwrap.dedent("""\
                        bedfile with transcript coordinates to convert
                        \n"""),
                        required = True)

    args = parser.parse_args() 
    
    bedgraph_fh = open(args.bedgraph)
    
    tbx = pysam.TabixFile(args.gene_bed)
    output_dir = "tmp"
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

    for tx in parse_bed(bedgraph_fh):
        print(tx)
        fout = open(os.path.join(output_dir, tx[3] + "_coverage.map"), 'w')
        for row in tbx.fetch(tx[0], tx[1], tx[2], parser=pysam.asBed()):
            print(row)
            new_start = max(row.start - tx[1], 1) 
            new_end = min(row.end - tx[1], tx[2] - tx[1] + 1)
            
            if new_start >= new_end:
                print("something is wrong with row " + str(row))
            
            aux_fields = [str(x) for x in row[3:]]
            fout.write("\t".join([tx[3], 
                                  str(new_start), 
                                  str(new_end),
                                  "\t".join(aux_fields)]) + "\n")
        fout.close()
 
if __name__ == '__main__': main()

