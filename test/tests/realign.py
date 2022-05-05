#   run as ctest .
#
# or
#
#   cd ../test ; env PYTHONPATH=../build python3 tests/realign.py ; cd ../build

import unittest
from pyvcflib import *

class RealignTest(unittest.TestCase):

    # Returns True or False.
    def test1(self):
        vcf = VariantCallFile()
        vcf.openFile("../samples/10158243.vcf")
        rec = Variant(vcf)
        vcf.getNextVariant(rec)
        self.assertEqual(rec.name,"grch38#chr4")
        self.assertEqual(rec.ref,'ACCCCCACCCCCACC')
        self.assertEqual(rec.alt,['ACC', 'AC', 'ACCCCCACCCCCAC', 'ACCCCCACC', 'ACA'])
        newvcf = rec.legacy_parsedAlternates(False,False,False,10.0,-9.0,15.0,6.66,0.0,"","",False,True)
        # for key, value in newvcf.items():
        #     print(f'Key: {key}: ')
        #     for a in value:
        #         print(f'{a.repr} ')

    def test_wftailbug(self):
        vcf = VariantCallFile()
        vcf.openFile("../samples/10134514.vcf")
        rec = Variant(vcf)
        vcf.getNextVariant(rec)
        self.assertEqual(rec.name,"grch38#chr4_10083863-10181258.vcf:grch38#chr4")
        self.assertEqual(rec.ref,'GGAGAATCCCAATTGATGG')
        self.assertEqual(rec.alt,['GTAGCATCCCAAGTGATGT', 'GTAGAATCCCAATTGATGT', 'GGAGCATCCCAATTGATGG', 'GG'])
        sw = rec.legacy_parsedAlternates(False,False,False,10.0,-9.0,15.0,6.66,0.0,"","",False,True)
        for key, value in sw.items():
            print(f'SW allele key: {key}: ')
            for a in value:
                print(f'SW allele: {a.repr} ')
        # note wf ignores paramaters
        wf = rec.legacy_parsedAlternates(False,False,False,10.0,-9.0,15.0,6.66,0.0,"","",True,True)
        for key, value in sw.items():
            print(f'WF allele key: {key}: ')
            for a in value:
                print(f'WF allele: {a.repr} ')


if __name__ == '__main__':
    unittest.main()
