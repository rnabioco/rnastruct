import gzip
import os
import mpileup_to_counts
import utils
import subprocess
import sys

class TestParsingPileup:
    "read pileup format from test file and test"


class TestPileuptocounts:

  def test_indels_pos(self, tmp_path):
      bam = "test_data/pileup_to_counts/chr16_fus.bam"
      fasta = "test_data/fasta/chr12_16_17_subset.fa.gz"
      out =  tmp_path / "fus_counts.tsv.gz"
      min_depth_arg = 1
      expected_output = "test_data/mpileup_to_counts/expected_pileup_to_counts_pos_overlaps.tsv"
      
      mpileup_args = utils.conv_args(utils.get_pileup_args())
      mpileup_args = mpileup_args.split()
      pileup_cmd =  ["samtools", "mpileup", "-f", fasta] + \
                  mpileup_args +  [bam]
      pileup_run = subprocess.Popen(pileup_cmd,
            shell=False,
            stderr = sys.stderr,
            stdout = subprocess.PIPE,
            universal_newlines = True) 
      
      fout = gzip.open(out, 'wt')
      mpileup_to_counts.mpileup_to_counts(pileup_run.stdout, 
            fout, 
            min_depth = min_depth_arg, 
            return_comp = False,
            debug = False)
      fout.close()      
      with open(expected_output, 'r') as exp:
          with gzip.open(out, 'rt') as res:
            for exp_idx, exp_line in enumerate(exp): 
              line_found = False
              for res_idx,res_line in enumerate(res):
                  if exp_line.split() != res_line.split():
                      continue
                  else:
                      line_found = True
                      break
              if not line_found:
                  msg = "line {}, not found in output\n{}".format(exp_idx,
                          " ".join(exp_line.split()))
                  assert False , msg 
              
  def test_indels_neg(self, tmp_path):
      bam = "test_data/pileup_to_counts/chr17_top2a.bam"
      fasta = "test_data/fasta/chr12_16_17_subset.fa.gz"
      out =  tmp_path / "top2a_counts.tsv.gz"
      min_depth_arg = 1
      expected_output = "test_data/mpileup_to_counts/expected_pileup_to_counts_neg_overlaps.tsv"
      
      mpileup_args = utils.conv_args(utils.get_pileup_args())
      mpileup_args = mpileup_args.split()
      pileup_cmd =  ["samtools", "mpileup", "-f", fasta] + \
                  mpileup_args +  [bam]
      pileup_run = subprocess.Popen(pileup_cmd,
            shell=False,
            stderr = sys.stderr,
            stdout = subprocess.PIPE,
            universal_newlines = True) 
      
      fout = gzip.open(out, 'wt')
      mpileup_to_counts.mpileup_to_counts(pileup_run.stdout, 
            fout,
            min_depth = min_depth_arg, 
            return_comp = True,
            debug = False)

      fout.close()      
      
      with open(expected_output, 'r') as exp:
          with gzip.open(out, 'rt') as res:
            for exp_idx, exp_line in enumerate(exp): 
              line_found = False
              for res_idx,res_line in enumerate(res):
                  if exp_line.split() != res_line.split():
                      continue
                  else:
                      line_found = True
                      break
              if not line_found:
                  msg = "line {}, not found in output\n{}".format(exp_idx,
                          " ".join(exp_line.split()))
                  assert False , msg 
