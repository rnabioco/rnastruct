import gzip
import os
import sys
import subprocess
import shutil
import pysam

class TestSimple:
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
      
      ivls_to_check = ["chr12", "6646357", "+", 103]
      ivls_found = False
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
            
            if [chrom,pos,strand,depth] == ivls_to_check:
                ivls_found = True
            

      strands = list(strands)
      chroms = list(chroms)
      assert idx > 1 
      assert len(strands) == 2 
      assert len(chroms) == 2 and all([x in chroms for x in ["chr12", "chr17"]]) 
      
      assert ("chr12", "+") in chr_strands
      assert ("chr17", "-") in chr_strands
      assert ivls_found

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

class TestComplex:
    bam = "test_data/chr16_fus_dms.bam"
    fa = "test_data/fasta/chr12_16_17_subset.fa.gz"

    def test_fus(self, tmp_path):
      outpre = tmp_path / "fus_"
      outpre = str(outpre)
      if os.path.exists("fus_test/"):
        shutil.rmtree("fus_test/", ignore_errors=False, onerror=None)
      exe_dir = os.path.dirname(os.path.abspath((__file__)))
      cmd_args = [
              os.path.join(exe_dir, "base_counter.py"),
              "-b", self.bam,
              "-f", self.fa,
              "-d", "10",
              "-o", outpre,
              "-k", "-v", 
              "-n", "ATCG"]
      out = subprocess.run(cmd_args, stdout = sys.stdout, stderr = sys.stderr)
      assert out.returncode == 0
      tf = pysam.TabixFile(outpre + "pileup_table.tsv.bgz")
      for rec in tf.fetch("chr16:31199224-31199224"):
        # indel is counted
        assert rec.split("\t")[11] == "1"
      
      for rec in tf.fetch("chr16:31199225-31199225"):
        # mismatches are counted
        assert rec.split("\t")[10] == "3"
      
      for rec in tf.fetch("chr16:31199225-31199225"):
        # mismatches are counted
        assert rec.split("\t")[10] == "3"


