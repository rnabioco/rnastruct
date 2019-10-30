import gzip
import os
import pileup_to_counts

class TestPileuptocounts:
    
  def test_fwd(self, tmp_path):
    vcf = "test_data/small_fwd.vcf"
    out = tmp_path / "counts.tsv.gz"
    min_depth_arg = 100
    
    pileup_to_counts.vcf_to_counts(vcf, 
            out, 
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
    out = tmp_path / "counts.tsv.gz"
    min_depth_arg = 100
    
    pileup_to_counts.vcf_to_counts(vcf, 
            out, 
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

