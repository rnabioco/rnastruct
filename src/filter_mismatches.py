#!/usr/bin/env python3
import argparse
import os, sys
import pysam

"""
Parse bam file and exclude alignments with non-A-G or T-C mismatches
"""


def count_bad_mismatches(alignment, good_mismatches):
    """ check for mismatches in get_aligned_pairs() output
    report number of non-good mismatches
    """
    
    bad_mismatches = 0
    
    seq = alignment.query_sequence

    try:
      aligned_pairs = alignment.get_aligned_pairs(with_seq = True,
              matches_only = False)
    except ValueError:
      sys.exit("MD flag not found in bamfile, use samtools calmd to add MD flag")
    
    for pos in aligned_pairs:
      if pos[0] is None or pos[2] is None:
          continue
      ref = pos[2]
      query = seq[pos[0]]
      
      if ref.islower():
        """ mismatches are encoded as lowercase ref nucleotide"""
        good_mismatch = False
        for mismatch_pair in good_mismatches:
          if ref == mismatch_pair[0].lower() and query == mismatch_pair[1]:
            good_mismatch = True
            break
        if good_mismatch:
            continue
        else:
            bad_mismatches += 1
    
    return bad_mismatches

def process_bam(bam, output_fn, mm_allowance, mm_proportion, mismatch_nts):
    
    bamfile_obj = pysam.AlignmentFile(bam, "rb" )

    outfile_obj = pysam.AlignmentFile(output_fn, 'wb', 
            template = bamfile_obj)

    for alignment in bamfile_obj.fetch():
        
        if alignment.is_unmapped:
            continue
        
        non_ag_mm = count_bad_mismatches(alignment, mismatch_nts)    
        
        if mm_proportion:
            mm_freq = float(non_ag_mm) / alignment.query_length
            if mm_freq > mm_allowance:
              continue
        else:
            if non_ag_mm > mm_allowance:
              continue
        
        outfile_obj.write(alignment)
        
    bamfile_obj.close()
    outfile_obj.close()

def main():
    
    parser = argparse.ArgumentParser(description="""
    Remove alignments with non A-G or T-C mismatches and return a bam
    file. The input bam file must have MD flags.""" )

    parser.add_argument('-b',
                          '--bam',
                          help ='indexed bam file input',
                       required = True)

    parser.add_argument('-o',
                          '--output_bam',
                          help =""" 
                          bam output file name
                          """,
                          required = True)
    parser.add_argument('-m',
                          '--mismatch_allowance',
                          help =""" Number of allowed non-A-G or T-C
                          mismatches, or proportion of read length allows
                          to have non-A-G or T-C mismatches. If argument
                          is >= 1 then it will be assumed to be the allowed
                          count of mismatches. If the argument is a
                          float between 0 and 1 then it will be assumed to
                          be a proportion. 
                          i.e 
                          -m 1 `allow alignments to have 1 non-A-G`
                          -m 0.01 `allow alignments to have 1%% non-A-G`
                          """,
                          required = False, 
                          default = 1,
                          type = float)
    parser.add_argument('-n',
                          '--mismatch_nts',
                          help =""" 
                          mismatch types allowed, defaults to A to G and T
                          to C. Supply each mismatch pair as strings with reference and
                          mismatch. Separate each mismatch pair with a
                          space.
                          i.e. "AG" "TC" will keep A to G and
                          T to C mismatches.
                          """,
                          required = False,
                          nargs = "+",
                          default = ["AG", "TC"])

    args=parser.parse_args()
    
    bam_name = args.bam
    output_bam = args.output_bam
    mm_count = args.mismatch_allowance
    nts = args.mismatch_nts

    if mm_count >= 1 or mm_count == 0:
      mm_proportion = False
      mm_allowances = int(mm_count)
    elif mm_count > 0 and mm_count < 1:
      mm_proportion = True
      mm_allowances = mm_count
    else:
      sys.exit("unable to parse supplied argument -m {}".format(mm_count))

    process_bam(bam_name, output_bam, mm_allowances, mm_proportion, nts) 
    
if __name__ == '__main__': main()
