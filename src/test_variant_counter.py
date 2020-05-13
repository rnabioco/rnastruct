import gzip
import os
import variant_counter
import pysam
import yaml
from cyvcf2 import VCF
from yaml import load, dump

class TestPosStrand:
    bam_fn = "test_data/variant_counter/chr16_fus_dms_delfixed.bam"
    vcf_fn = "test_data/variant_counter/fus.vcf.gz"
    config_fn = "pileup.yaml"

    vcf = VCF(vcf_fn)
    bam = pysam.AlignmentFile(bam_fn, 'rb')
  
    with open(config_fn) as f:
        config_data = f.read()

    pup_args = load(config_data)
    
    #def test_basic_usage(self):
    #    bam_name = variant_counter.scan_bam(self.bam_fn, self.vcf_fn, True, self.pup_args)
    
    def test_count_variants(self):
        
        for v in self.vcf("chr16:31196456-31196457"):
            """
            basic usage
            """
            if v.INFO.get('INDEL') is None:
                continue
            pileup_itr = self.bam.pileup(v.CHROM, v.start,
                    v.start + 1,
                    truncate = True,
                    **self.pup_args)
            res = variant_counter.count_variants(pileup_itr, v)
            
            assert len(res.keys()) == 1
            assert list(res.keys())[0] == "AC"
            assert res["AC"] == 2
        
        for v in self.vcf("chr16:31199093-31199093"):
            """
            multiple indels
            """
            if v.INFO.get('INDEL') is None:
                continue
            pileup_itr = self.bam.pileup(v.CHROM, v.start,
                    v.start + 1,
                    truncate = True,
                    **self.pup_args)
            res = variant_counter.count_variants(pileup_itr, v)
            
            assert len(res.keys()) == 2 
            assert all(x in res.keys() for x in ["TAAAA", "TAAA"])
            assert res["TAAAA"] == 3
            assert res["TAAA"] == 1
        
        for v in self.vcf("chr16:31199185-31199185"):
            """
            new allele discovered and multiple indels
            """
            if v.INFO.get('INDEL') is None:
                continue
            pileup_itr = self.bam.pileup(v.CHROM, v.start,
                    v.start + 1,
                    truncate = True,
                    **self.pup_args)
            res = variant_counter.count_variants(pileup_itr, v)
            assert len(res.keys()) == 4 
            assert all(x in res.keys() for x in ["TAA", "TAAAAAAA", "TAAA"])
            assert "T-1N" in res.keys()
            assert res["T-1N"] == 10
            assert res["TAAAAAAA"] == 4
            assert res["TAAA"] == 2
            assert res["TAA"] == 2
