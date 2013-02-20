"""
Some simple unit tests for the vcfnp extension.

"""


from vcfnp import variants
from nose.tools import eq_


def test_variants():
    a = variants('sample.vcf')
    eq_(9, len(a))
    eq_('19', a[0]['CHROM'])
    eq_(111, a[0]['POS'])
    eq_('rs6054257', a[2]['ID'])
    eq_('A', a[0]['REF'])
    eq_('ATG', a[8]['ALT'][1])
    eq_(10.0, a[1]['QUAL'])
    eq_(True, a[2]['FILTER']['PASS'])
    eq_(2, a[0]['num_alleles'])
    eq_(False, a[5]['is_snp'])
    
#array([ ('19', 111, '.', 'A', ['C', ''], 9.600000381469727, (False, False, False), 2, True),
#       ('19', 112, '.', 'A', ['G', ''], 10.0, (False, False, False), 2, True),
#       ('20', 14370, 'rs6054257', 'G', ['A', ''], 29.0, (True, False, False), 2, True),
#       ('20', 17330, '.', 'T', ['A', ''], 3.0, (False, True, False), 2, True),
#       ('20', 1110696, 'rs6040355', 'A', ['G', 'T'], 67.0, (True, False, False), 3, True),
#       ('20', 1230237, '.', 'T', ['.', ''], 47.0, (True, False, False), 2, False),
#       ('20', 1234567, 'microsat1', 'G', ['GA', 'GAC'], 50.0, (True, False, False), 3, False),
#       ('20', 1235237, '.', 'T', ['.', ''], 0.0, (False, False, False), 2, False),
#       ('X', 10, 'rsTest', 'AC', ['A', 'ATG'], 10.0, (True, False, False), 3, False)], 
#      dtype=[('CHROM', '|S12'), ('POS', '<i4'), ('ID', '|S12'), ('REF', '|S12'), ('ALT', '|S12', (2,)), ('QUAL', '<f4'), ('FILTER', [('PASS', '|b1'), ('q10', '|b1'), ('s50', '|b1')]), ('num_alleles', '|u1'), ('is_snp', '|b1')])


def test_variants_region():
    a = variants('sample.vcf.gz', region='20')
    eq_(6, len(a))
    
    
def test_variants_count():
    a = variants('sample.vcf', count=3)
    eq_(3, len(a))

    
    
