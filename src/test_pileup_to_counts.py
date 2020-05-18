import gzip
import os
import pileup_to_counts
import utils

class TestPileuptocounts:
  
  def test_fwd(self, tmp_path):
    vcf = "test_data/small_fwd.vcf"
    bam = "test_data/small.fwd.bam"
    out = tmp_path / "counts.tsv.gz"
    min_depth_arg = 100
    
    pileup_to_counts.vcf_to_counts(vcf, 
            bam,
            out, 
            utils.get_pileup_args(),
            min_depth = min_depth_arg, 
            return_comp = False,
            debug = False)
  
    assert os.path.isfile(out)
    
    n_lines = 0
    strand = set()
    chroms = set()
    n_fields = 0
    min_depth = min_depth_arg
    cols = []
  
    with gzip.open(out, 'rt') as f:
      for idx,line in enumerate(f):
          
          n_lines += 1
          
          fields = line.strip().split("\t")
          if idx == 0:
            n_fields = len(fields)
            cols = fields
            continue
  
          chroms.add(fields[0])
          strand.add(fields[2])
          
          fields[4] = int(fields[4])
          if fields[4] < min_depth:
              min_depth = fields[4]

    strand = list(strand)
    chroms = list(chroms)
    assert n_lines > 1 
    assert len(strand) == 1 and strand[0] == "+" 
    assert len(chroms) == 1 and all([x in chroms for x in ["chr12"]]) 
    assert min_depth >= min_depth_arg

  def test_rev(self, tmp_path):
    vcf = "test_data/small_rev.vcf"
    bam = "test_data/small.rev.bam"
    out = tmp_path / "counts.tsv.gz"
    min_depth_arg = 100
    
    pileup_to_counts.vcf_to_counts(vcf, 
            bam,
            out, 
            utils.get_pileup_args(),
            min_depth = min_depth_arg, 
            return_comp = True,
            debug = False)
  
    assert os.path.isfile(out)
    
    n_lines = 0
    strand = set()
    chroms = set()
    n_fields = 0
    min_depth = min_depth_arg
    cols = []
  
    with gzip.open(out, 'rt') as f:
      for idx,line in enumerate(f):
          
          n_lines += 1
          
          fields = line.strip().split("\t")
          if idx == 0:
            n_fields = len(fields)
            cols = fields
            continue
  
          chroms.add(fields[0])
          strand.add(fields[2])
          
          fields[4] = int(fields[4])
          if fields[4] < min_depth:
              min_depth = fields[4]

    strand = list(strand)
    chroms = list(chroms)
    assert n_lines > 1 
    assert len(strand) == 1 and strand[0] == "-" 
    assert len(chroms) == 1 and all([x in chroms for x in ["chr17"]]) 
    assert min_depth >= min_depth_arg
  

  def test_indels_pos(self, tmp_path):
      vcf = "test_data/pileup_to_counts/fus.overlaps.vcf.gz"
      bam = "test_data/pileup_to_counts/chr16_fus.bam"
      out =  tmp_path / "fus_counts.tsv.gz"
      min_depth_arg = 1
      expected_output = "test_data/pileup_to_counts/expected_pileup_to_counts_pos_overlaps.tsv"
      
      pileup_to_counts.vcf_to_counts(vcf, 
            bam,
            out, 
            utils.get_pileup_args(custom_args = {"ignore_overlaps":False}),
            min_depth = min_depth_arg, 
            return_comp = False,
            debug = False)
      
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
      vcf = "test_data/pileup_to_counts/top2a.overlaps.vcf.gz"
      bam = "test_data/pileup_to_counts/chr17_top2a.bam"
      out =  tmp_path / "top2a_counts.tsv.gz"
      min_depth_arg = 10
      expected_output = "test_data/pileup_to_counts/expected_pileup_to_counts_neg_overlaps.tsv"
      
      pileup_to_counts.vcf_to_counts(vcf, 
            bam,
            out, 
            utils.get_pileup_args(custom_args = {"ignore_overlaps":False}),
            min_depth = min_depth_arg, 
            return_comp = True,
            debug = False)
      
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
              
