import pytest
from utils import is_complex_indel 

class TestComplexVariants:

    def test_complex_variants(self):

        assert not is_complex_indel("AAT", "AT")
        assert not is_complex_indel("TATATA", "TA")
        assert not is_complex_indel("GC", "G")
        assert not is_complex_indel("GGG", "GGGG")
        assert is_complex_indel("AAT", "ATCGCGAGT")
        assert is_complex_indel("CCGCGT", "CCCA")
    
    def test_MVP_throws_error(self):
        with pytest.raises(Exception):
            is_complex_indel("GGAAGG", "GTCAAT")
    
    def test_snp_throws_error(self):
        with pytest.raises(Exception):
            is_complex_indel("G", "C")

