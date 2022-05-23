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
        self.assertEqual(info['AC'],['11', '7', '1', '3'])
        unique = {}
        a = None
        for alt0, value in wf.items():
            is_rev = value[1]
            for matches in value[0]:
                ref = matches.ref
                aligned = matches.alt
                tag = f'{alt0}:{matches.position}:{ref}/{aligned}'
                if rec.ref == aligned:
                    alt_index = -1
                    AC = None
                    AF = None
                    AN = None
                else:
                    alt_index = rec.alt.index(alt0) # Raises a ValueError if there is no such item
                    AC = int(info['AC'][alt_index])
                    AF = float(info['AF'][alt_index])
                    AN = int(info['AN'][0])
                relpos = matches.position - rec.pos
                unique[tag] = {
                    'pos0': rec.pos,
                    'ref0': rec.ref,
                    'alt0': alt0,
                    'ref1': ref,
                    'algn': aligned,
                    'pos1': matches.position,
                    'altidx': alt_index, # zero based
                    'relpos': relpos,
                    'AC': AC,
                    'AF': AF,
                    'AN': AN,
                    'is_rev': is_rev}
        # Did we get all?
        self.assertEqual(len(unique.items()),18)
        # Display
        uniqsorted = sorted(unique.items(),key = lambda r: r[1]['pos1'])
        print(json.dumps(uniqsorted,indent=4))
        # Check if all alleles were used by counting 'altidx'
        idxs = set(map(lambda k: unique[k]['altidx'],unique.keys()))
        self.assertEqual(len(idxs),len(rec.alt)+1)
        # We are now going to merge records using a dict. From
        #
        # 10134532:AGAATCCCAATTGATGG/
        # 10134532:A/C
        # 10134532:G/T
        # 10134532:G/T
        # 10134532:A/C
        # 10134532:T/G
        variants = {}
        for k,v in uniqsorted:
            ref = v['ref1']
            aligned = v['algn']
            if ref != aligned:
                ntag = f"{v['pos1']}:{ref}/{aligned}"
                print(f"{ntag} AC={v['AC']}")
                if ntag in variants:
                    variants[ntag]['AC'] += v['AC']
                    # Check AN number is equal so we can compute AF by addition
                    self.assertEqual(variants[ntag]['AN'],v['AN'])
                    variants[ntag]['AF'] += v['AF']
                else:
                    variants[ntag] = v
        # print(variants)
        # print(json.dumps(variants,indent=4))
        self.assertEqual(variants['10134532:G/T']['AC'],18)
        for key in variants:
            v = variants[key]
            ref_len = len(v['ref1'])
            aln_len = len(v['algn'])
            type = None
            if aln_len < ref_len:
                type = 'del'
            elif aln_len > ref_len:
                type = 'del'
            elif aln_len == ref_len:
                if ref_len == 1:
                    type = 'snp'
                else:
                    type = 'mnp'
            variants[key]['type'] = type
            # Set origin
            variants[key]['origin'] = f"{rec.name}:{rec.pos}"
        print(json.dumps(variants,indent=4))
        samples = rec.samples
        gts = []
        for name in rec.sampleNames:
            # print(name,samples[name])
            gt = (samples[name]['GT'])[0].split("|")
            # print(gt)
            gts.append(list(map(lambda item: int(item) if item.isdigit() else None,gt)))
        print(gts)
        # for each variant translate genotypes
        for key in variants:
            print(key)
            rec = variants[key]
            idx1 = rec['altidx']+1
            print(idx1)
            genotypes = []
            for gt in gts:
                # print(gt)
                genotypes.append(list(map(lambda item: item,gt)))
            print(list(genotypes))
            # Now we neet to plug in the new indices
            for gt in genotypes:
                if gt[0] == idx1:
                    gt[0] = 1
                else:
                    if gt[0] != None:
                        gt[0] = 0
                if len(gt)>1:
                    if gt[1] == idx1:
                        gt[1] = 1
                    else:
                        if gt[1] != None:
                            gt[1] = 0
            print("---")
            print(list(genotypes))
# [[0, 1], [1, 0], [0, 0], [0, 1], [0, 0], [1, 0], [1, 0], [1, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, None], [0, 0], [2, 2], [0, 0], [4, 0], [0, 0], [0, 1], [0, 1], [0, 1], [0, 2], [0, 0], [4, 0], [0, 2], [0, 0], [0, 0], [2, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [2, 0], [4, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 3], [0, 0], [0, 2], [0, 0], [1]]

# Note we are '2'

#    0|1     1|0     0|0     0|1     0|0     1|0     1|0     1|0     0|0    0|0     0|0     0|0     0|0     0|.          0|0    1|1      0|0    .|0     0|0     0|1     0|1     0|1     0|1     0|0     .|0     0|1     0|0     0|0  1|0     0|0     0|0     0|0     0|0     0|0     1|0     .|1     0|0     0|0     0|0     0|0     0|0     0|0  0|1     0|0     1


if __name__ == '__main__':
    unittest.main()
