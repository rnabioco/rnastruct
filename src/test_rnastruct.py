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


class TestBaseCounter:
    """ testing output file content for now """
    fwd_pre = "test_data/base_counter/fwd_"
    rev_pre = "test_data/base_counter/rev_"
    both_pre = "test_data/base_counter/both_"

    def test_both(self):
      
      n_lines = 0
      strands = set()
      chroms = set()
      n_fields = 0
      min_depth = 100 
      cols = []
      chr_strands = set()
      uid_set = set()

      with gzip.open(self.both_pre + "pileup_table.tsv.bgz", 'rt') as f:
        for idx,line in enumerate(f):
            fields = line.strip().split("\t")
            if idx == 0:
              n_fields = len(fields)
              cols = fields
              continue
  
            chrom = fields[0]
            pos = fields[1]
            strand = fields[2]
            chroms.add(chrom)
            strands.add(strand)
            chr_strands.add((chrom, strand))
       
            # see if any positions are duplicated
            if (chrom, pos, strand) in uid_set:
                raise AssertionError("chrom, pos, and strand shouldn't be duplicated")
            else:
                uid_set.add((chrom, pos, strand))

            depth = int(fields[4])
            if depth < min_depth:
                raise AttributeError("depth of {} detected".format(depth))

      strands = list(strands)
      chroms = list(chroms)
      assert idx > 1 
      assert len(strands) == 2 
      assert len(chroms) == 2 and all([x in chroms for x in ["chr12", "chr17"]]) 
      
      assert ("chr12", "+") in chr_strands
      assert ("chr17", "-") in chr_strands

    def test_fwd(self):
        strands = set()
        with gzip.open(self.fwd_pre + "pileup_table.tsv.bgz", 'rt') as f:
          for idx,line in enumerate(f):
            fields = line.strip().split("\t")
            if idx == 0:
              continue  
            strands.add(fields[2])

          assert len(strands) == 1
    
    def test_rev(self):
        strands = set()
        with gzip.open(self.rev_pre + "pileup_table.tsv.bgz", 'rt') as f:
          for idx,line in enumerate(f):
            fields = line.strip().split("\t")
            if idx == 0:
              continue  
            strands.add(fields[2])

          assert len(strands) == 1
