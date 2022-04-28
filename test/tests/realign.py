# env PYTHONPATH=../build python3 tests/realign.py

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
        for key, value in newvcf.items():
            print(f'Key: {key}: ')
            for a in value:
                print(f'{a.repr} ')


if __name__ == '__main__':
    unittest.main()
