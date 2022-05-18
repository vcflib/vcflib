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
        # A dict is returned of alleles with variants and is_reversed
        wf = rec.parsedAlternates(False,True,False,"","",wfa_params,True,64,True)
        print(f'ref={rec.ref}')
        print(rec.info)
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
        info = rec.info
        unique = {}
        for alt0, value in wf.items():
            is_rev = value[1]
            for matches in value[0]:
                ref = matches.ref
                aligned = matches.alt
                tag = f'{alt0}:{a.position}:{ref}/{aligned}'
                if rec.ref == aligned:
                    alt_index = -1
                    AC = None
                    AF = None
                else:
                    alt_index = rec.alt.index(alt0) # Raises a ValueError if there is no such item
                    AC = info['AC'][alt_index]
                    AF = info['AF'][alt_index]
                relpos = a.position - rec.pos
                unique[tag] = {
                    'pos0': rec.pos,
                    'ref0': rec.ref,
                    'alt0': alt0,
                    'ref1': ref,
                    'algn': aligned,
                    'pos': a.position, # points where exactly? FIXME
                    'altidx': alt_index,
                    'relpos': relpos,
                    'AC': AC,
                    'AF': AF,
                    'is_rev': is_rev}
        # uniqsorted = sorted(unique.items(),key = lambda r: r[1]['pos'])
        uniqsorted = unique
        print(json.dumps(uniqsorted,indent=4))
        self.assertEqual(len(uniqsorted),16)
        # Check if all alleles were used by counting 'altidx'
        idxs = set(map(lambda k: unique[k]['altidx'],unique.keys()))
        self.assertEqual(len(idxs),len(rec.alt)+1)
        # We are now going to adjust the info fields
        self.assertEqual(rec.info['AC'],['11', '7', '1', '3'])

if __name__ == '__main__':
    unittest.main()
