#   run as ctest .
#
# or
#
#   cd ../test ; env PYTHONPATH=../build python3 tests/realign.py ; cd ../build

import unittest
from pyvcflib import *
import json

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
        sw = rec.legacy_parsedAlternates(False,False,False,10.0,-9.0,15.0,6.66,0.0,"","",False,False)

    def test_sw_wf_compare(self):
        vcf = VariantCallFile()
        vcf.openFile("../samples/10134514.vcf")
        rec = Variant(vcf)
        vcf.getNextVariant(rec)
        self.assertEqual(rec.name,"grch38#chr4_10083863-10181258.vcf:grch38#chr4")
        self.assertEqual(rec.ref,'GGAGAATCCCAATTGATGG')
        self.assertEqual(rec.alt,['GTAGCATCCCAAGTGATGT', 'GTAGAATCCCAATTGATGT', 'GGAGCATCCCAATTGATGG', 'GG'])
        sw = rec.legacy_parsedAlternates(False,False,False,10.0,-9.0,15.0,6.66,0.0,"","",False,False)
        # for key, value in sw.items():
        #     print(f'SW allele key: {key}: ')
        #     for a in value:
        #         print(f'               {a.position}/{a.ref}/{a.alt} ')
        # note wf ignores paramaters
        wf = rec.legacy_parsedAlternates(False,False,False,10.0,-9.0,15.0,6.66,0.0,"","",True,False)
        # for key, value in wf.items():
        #     print(f'WF allele key: {key}: ')
        #     for a in value:
        #        print(f'               {a.position}/{a.ref}/{a.alt} ')
        self.assertEqual(len(wf),5)
        self.assertEqual(len(wf),len(sw))
        self.assertEqual(sw['GGAGAATCCCAATTGATGG'][0].alt,wf['GGAGAATCCCAATTGATGG'][0].alt)

    def test_wfbug2(self):
        vcf = VariantCallFile()
        vcf.openFile("../samples/10134514.vcf")
        rec = Variant(vcf)
        vcf.getNextVariant(rec)
        self.assertEqual(rec.name,"grch38#chr4_10083863-10181258.vcf:grch38#chr4")
        self.assertEqual(rec.ref,'GGAGAATCCCAATTGATGG')
        self.assertEqual(rec.alt,['GTAGCATCCCAAGTGATGT', 'GTAGAATCCCAATTGATGT', 'GGAGCATCCCAATTGATGG', 'GG'])
        wfa_params = wavefront_aligner_attr_default
        # string paramString = "0,19,39,3,81,1";
        wfa_params.distance_metric = distance_meric_t.gap_affine_2p
        wfa_params.affine2p_penalties.match = 0
        wfa_params.affine2p_penalties.mismatch = 19
        wfa_params.affine2p_penalties.gap_opening1 = 39
        wfa_params.affine2p_penalties.gap_extension1 = 3
        wfa_params.affine2p_penalties.gap_opening2 = 81
        wfa_params.affine2p_penalties.gap_extension2 = 1
        wfa_params.alignment_scope = alignment_scope_t.compute_alignment;
        wf = rec.parsedAlternates(False,True,False,"","",wfa_params,True,64,True)
        for key, value in wf.items():
            print(f'WF2 allele key: {key}: ')
            for a in value[0]:
                print(f'               {a.position}:{a.ref}/{a.alt} ')
        self.assertEqual(len(wf),5)
        gg0 = wf['GG'][0][0]
        gg1 = wf['GG'][0][1]
        self.assertEqual(gg0.alt,"GG")
        self.assertEqual(gg1.alt,"")
        self.assertEqual(wf['GGAGAATCCCAATTGATGG'][0][0].alt,"GGAGAATCCCAATTGATGG")
        # collect unique alleles
        unique = {}
        alt_index = -1
        for alt1, value in wf.items():
            alt_index += 1
            for a in value[0]:
                ref = a.ref
                alt = a.alt
                if ref != alt and alt != "":
                    tag = f'{alt1}:{a.position}:{ref}/{alt}'
                    unique[tag] = { 'alt1': alt1,
                                    'pos': a.position,
                                    'ref': ref,
                                    'aidx': alt_index,
                                    'alt': alt}
        print(json.dumps(unique,indent=4))


if __name__ == '__main__':
    unittest.main()
