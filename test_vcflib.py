"""
Some simple unit tests for the Cython vcflib wrapper.

"""


from vcflib import PyVariantCallFile, TYPE_INTEGER, TYPE_FLOAT, TYPE_STRING, TYPE_BOOL
from nose.tools import eq_


def test_len():
    vcf = PyVariantCallFile('sample.vcf')
    eq_(9, len(vcf))
    

def test_metadata():
    vcf = PyVariantCallFile('sample.vcf')
    eq_(['NS', 'AN', 'AC', 'DP', 'AF', 'AA', 'DB', 'H2'], vcf.infoIds)
    eq_(['GT', 'GQ', 'DP', 'HQ'], vcf.formatIds)
    eq_(['q10', 's50'], vcf.filterIds)
    eq_(TYPE_INTEGER, vcf.infoTypes['NS'])
    eq_(TYPE_FLOAT, vcf.infoTypes['AF'])
    eq_(TYPE_STRING, vcf.infoTypes['AA'])
    eq_(TYPE_BOOL, vcf.infoTypes['DB'])
    eq_(TYPE_STRING, vcf.formatTypes['GT'])
    eq_(TYPE_INTEGER, vcf.formatTypes['GQ'])
    eq_(1, vcf.infoCounts['AA'])
    eq_(2, vcf.formatCounts['HQ'])
    eq_(['NA00001', 'NA00002', 'NA00003'], vcf.sampleNames)


def test_fixed_fields():
    vcf = PyVariantCallFile('sample.vcf')
    v = iter(vcf).next() # first variant
    eq_('19', v.CHROM)
    eq_(111, v.POS)
    eq_('.', v.ID)
    eq_('A', v.REF)
    eq_(['C'], v.ALT)
    eq_(9.6, v.QUAL)
    eq_(['.'], v.FILTER) # split in C++
    
    
def test_info():
    vcf = PyVariantCallFile('sample.vcf')
    v = list(iter(vcf))[4] # fifth variant
    expect = {'AA': ['T'], 'NS': ['2'], 'DP': ['10'], 'AF': ['0.333', '0.667']}
    eq_(expect, v.INFO)


def test_samples():
    vcf = PyVariantCallFile('sample.vcf')
    v = list(iter(vcf))[2] # third variant
    expect = {'GT': ['0|0'], 'HQ': ['51', '51'], 'GQ': ['48'], 'DP': ['1']}
    eq_(expect, v.samples['NA00001'])


def test_region():
    vcf = PyVariantCallFile('sample.vcf.gz')
    vcf.setRegion('20')
    eq_(6, len(vcf))
    vcf.setRegion('20:1000000-2000000')
    eq_(4, len(vcf))
    vcf.setRegion('20', 1000000, 2000000)
    eq_(4, len(vcf))
    
    