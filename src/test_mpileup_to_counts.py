import gzip
import os
import mpileup_to_counts
import utils
import subprocess
import sys
import pandas as pd

class TestParsingPileup:
    "read pileup format from test file and test"
    plp_fn = "test_data/mpileup_to_counts/pileups_to_test.txt"

    def test_basic(self, tmp_path):
        out_fn = tmp_path / "out.tsv.gz"
        out_fn = str(out_fn)
  
        f = open(self.plp_fn, 'r')
        fo = gzip.open(out_fn, 'wt')

        mpileup_to_counts.mpileup_to_counts(f,
                    fo,
                    min_depth = 10, 
                    max_del = 100,
                    header = True)
        f.close()
        fo.close()
        assert os.path.isfile(out_fn)
        
        df = pd.read_csv(out_fn, sep = "\t")
        
        assert all(df["strand"] == "+")
        assert all(df["depth"] > 10)
        assert all(df[df["pos"] == 31199145]["depth"] == 63)
        assert all(df[df["pos"] == 31199145]["refcount"] == 61)
        assert all(df[df["pos"] == 31199145]["tcount"] == 1)
        assert all(df[df["pos"] == 31199145]["mmcount"] == 1)
        # 7nt deletion reported
        assert all(df[df["pos"] == 31199146]["indelcount"] == 1)
        assert all(df[df["pos"] == 31199198]["indelcount"] == 17)
        assert all(df[df["pos"] == 31199197]["depth"] == 109)
        assert all(df[df["pos"] == 31199198]["depth"] == 106)
               
        assert all(df[df["pos"] == 31199185]["indelcount"] == 4)
        assert all(df[df["pos"] == 31199186]["indelcount"] == 10)
        assert all(df[df["pos"] == 31199187]["indelcount"] == 2)
        assert all(df[df["pos"] == 31199188]["indelcount"] == 2)

    def test_drop_deletions(self, tmp_path):
        out_fn = tmp_path / "out.tsv.gz"
        out_fn = str(out_fn)
  
        f = open(self.plp_fn, 'r')
        fo = gzip.open(out_fn, 'wt')

        mpileup_to_counts.mpileup_to_counts(f,
                    fo,
                    min_depth = 10, 
                    max_del = 4,
                    header = True)
        f.close()
        fo.close()
        assert os.path.isfile(out_fn)
        
        df = pd.read_csv(out_fn, sep = "\t")
        
        assert all(df["strand"] == "+")
        assert all(df["depth"] > 10)
        assert all(df[df["pos"] == 31199145]["depth"] == 62)
        assert all(df[df["pos"] == 31199145]["refcount"] == 61)
        assert all(df[df["pos"] == 31199145]["tcount"] == 1)
        assert all(df[df["pos"] == 31199145]["mmcount"] == 1)
        # 7nt not deletion reported
        assert all(df[df["pos"] == 31199146]["indelcount"] == 0)
        assert all(df[df["pos"] == 31199198]["indelcount"] == 17)
        assert all(df[df["pos"] == 31199197]["depth"] == 109)
        assert all(df[df["pos"] == 31199198]["depth"] == 106)


    def test_complement(self, tmp_path):
        out_fn = tmp_path / "out.tsv.gz"
        out_fn = str(out_fn)
  
        f = open(self.plp_fn, 'r')
        fo = gzip.open(out_fn, 'wt')

        mpileup_to_counts.mpileup_to_counts(f,
                    fo,
                    min_depth = 10,
                    return_comp = True,
                    max_del = 4,
                    header = True)
        f.close()
        fo.close()
        assert os.path.isfile(out_fn)
        
        df = pd.read_csv(out_fn, sep = "\t")
        
        assert all(df["strand"] == "-")
        assert all(df["depth"] > 10)
        assert all(df[df["pos"] == 31199145]["ref_base"] == "C")
        assert all(df[df["pos"] == 31199145]["depth"] == 62)
        assert all(df[df["pos"] == 31199145]["refcount"] == 61)
        assert all(df[df["pos"] == 31199145]["ccount"] == 61)
        assert all(df[df["pos"] == 31199145]["acount"] == 1)
        assert all(df[df["pos"] == 31199145]["mmcount"] == 1)
        # 7nt not deletion reported
        assert all(df[df["pos"] == 31199146]["indelcount"] == 0)
        assert all(df[df["pos"] == 31199197]["indelcount"] == 15)
        assert all(df[df["pos"] == 31199198]["indelcount"] == 2)
        assert all(df[df["pos"] == 31199197]["depth"] == 109)
        assert all(df[df["pos"] == 31199198]["depth"] == 106)
        
        assert all(df[df["pos"] == 31199185]["indelcount"] == 0)
        assert all(df[df["pos"] == 31199186]["indelcount"] == 18)
        assert all(df[df["pos"] == 31199187]["indelcount"] == 0)
        assert all(df[df["pos"] == 31199188]["indelcount"] == 0)

    def test_basic(self, tmp_path):
        out_fn = tmp_path / "out.tsv.gz"
        out_fn = str(out_fn)
  
        f = open(self.plp_fn, 'r')
        fo = gzip.open(out_fn, 'wt')

        mpileup_to_counts.mpileup_to_counts(f,
                    fo,
                    min_depth = 200, 
                    max_del = 4,
                    header = True)
        f.close()
        fo.close()
        assert os.path.isfile(out_fn)
        
        df = pd.read_csv(out_fn, sep = "\t")
        
        assert all(df["strand"] == "+")
        assert all(df["depth"] > 200)
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
