#   run as ctest .
#
# or
#
#   cd ../test ; env PYTHONPATH=../build python3 tests/realign.py ; cd ../build
#   cd ../test ; env PYTHONPATH=../build pytest tests/realign.py ; cd ../build

import unittest
from pyvcflib import *
import json

class RealignTest(unittest.TestCase):

    # Returns True or False.
    def test1(self):
        vcf = VariantCallFile()
        vcf.openFile("../samples/10158243.vcf")
        var = Variant(vcf)
        vcf.getNextVariant(var)
        self.assertEqual(var.name,"grch38#chr4")
        self.assertEqual(var.ref,'ACCCCCACCCCCACC')
        self.assertEqual(var.alt,['ACC', 'AC', 'ACCCCCACCCCCAC', 'ACCCCCACC', 'ACA'])
        sw = var.legacy_parsedAlternates(False,False,False,10.0,-9.0,15.0,6.66,0.0,"","",False,False)

    def test_sw_wf_compare(self):
        vcf = VariantCallFile()
        vcf.openFile("../samples/10134514.vcf")
        var = Variant(vcf)
        vcf.getNextVariant(var)
        self.assertEqual(var.name,"grch38#chr4_10083863-10181258.vcf:grch38#chr4")
        self.assertEqual(var.ref,'GGAGAATCCCAATTGATGG')
        self.assertEqual(var.alt,['GTAGCATCCCAAGTGATGT', 'GTAGAATCCCAATTGATGT', 'GGAGCATCCCAATTGATGG', 'GG'])
        sw = var.legacy_parsedAlternates(False,False,False,10.0,-9.0,15.0,6.66,0.0,"","",False,False)
        # for key, value in sw.items():
        #     print(f'SW allele key: {key}: ')
        #     for a in value:
        #         print(f'               {a.position}/{a.ref}/{a.alt} ')
        # note wf ignores paramaters
        wf = var.legacy_parsedAlternates(False,False,False,10.0,-9.0,15.0,6.66,0.0,"","",True,False)
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
        var = WfaVariant(vcf)
        vcf.getNextVariant(var)
        self.assertEqual(var.name,"grch38#chr4_10083863-10181258.vcf:grch38#chr4")
        self.assertEqual(var.ref,'GGAGAATCCCAATTGATGG')
        self.assertEqual(var.alt,['GTAGCATCCCAAGTGATGT', 'GTAGAATCCCAATTGATGT', 'GGAGCATCCCAATTGATGG', 'GG'])
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
        wf = var.wfa_parsedAlternates(True,True,False,"","",wfa_params,True,64,1,True)
        print(f'ref={var.ref}')
        print(var.info)
        for key1, value1 in wf.items():
            print(f'WF2 allele key: {key1}: ')
            for a in value1[0]:
                print(f'               {a.position}:{a.ref}/{a.alt} ')
        # Run a few tests
        self.assertEqual(len(wf),5)
        gg0 = wf['GG'][0][0]
        gg1 = wf['GG'][0][1]
        self.assertEqual(gg0.alt,"GG")
        self.assertEqual(gg1.alt,"G")
        self.assertEqual(wf['GGAGAATCCCAATTGATGG'][0][0].alt,"GGAGAATCCCAATTGATGG")
        # Collect unique alleles
        info = var.info
        self.assertEqual(info['AC'],['11', '7', '1', '3'])
        unique = {}
        a = None
        for alt0, wfvalue in wf.items(): # wfvalue is a compound of bool is_rev and alleles
            is_rev = wfvalue[1]
            for wfmatch in wfvalue[0]:
                ref = wfmatch.ref
                aligned = wfmatch.alt
                wfpos = wfmatch.position
                wftag = f'{alt0}:{wfpos}:{ref}/{aligned}'
                if var.ref == aligned:
                    alt_index = -1
                    AC = None
                    AF = None
                    AN = None
                else:
                    alt_index = var.alt.index(alt0) # Raises a ValueError if there is no such item
                    AC = int(info['AC'][alt_index])
                    AF = float(info['AF'][alt_index])
                    AN = int(info['AN'][0])
                relpos = wfpos - var.pos
                unique[wftag] = {
                    'pos0': var.pos,
                    'ref0': var.ref,
                    'alt0': alt0,
                    'ref1': ref,
                    'algn': aligned,
                    'pos1': wfpos,
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
        self.assertEqual(len(idxs),len(var.alt)+1)

        # Collect sample genotypes
        samples = var.samples
        gts = []
        for sname in var.sampleNames:
            # print(name,samples[name])
            gt = (samples[sname]['GT'])[0].split("|")
            # print(gt)
            gts.append(list(map(lambda item: int(item) if item.isdigit() else None,gt)))
        print(gts)
        # for each variant translate genotypes
        for tag,aln in uniqsorted:
            # print(tag)
            idx1 = aln['altidx']+1
            # print(idx1)
            genotypes = []
            for gt in gts:
                # print(gt)
                genotypes.append(list(map(lambda item: item,gt)))
            # print(list(genotypes))
            # Now we neet to plug in the new indices
            for gt in genotypes:
                for i,g in enumerate(gt):
                    if g == idx1:
                        gt[i] = 1 # only one genotype in play
                    else:
                        if g != None:
                            gt[i] = 0
            # print(list(genotypes))
            aln['samples'] = genotypes
        gts = None
        print(uniqsorted[10])
        self.assertEqual(uniqsorted[10][1]['samples'],[[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, None], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 1], [0, 0], [0, 0], [0, 0], [0]])

        # We are now going to merge records using a dict. From
        #
        # 10134515:G/T AC=7
        # 10134515:G/T AC=11
        # 10134516:AGAATCCCAATTGATGG/ AC=3
        # 10134518:A/C AC=1
        # 10134518:A/C AC=11
        # 10134526:T/G AC=11
        # 10134532:G/T AC=7
        # 10134532:G/T AC=11
        #   into
        # 10134515:G/T AC=18
        # 10134516:AGAATCCCAATTGATGG/ AC=3
        # 10134518:A/C AC=12
        # 10134526:T/G AC=11
        # 10134532:G/T AC=18
        variants = {} # store new hash
        k = None
        v = None
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
                    # Merge genotypes if they come from different alleles
                    if v['altidx'] != variants[ntag]['altidx']:
                        for i,samplesi in enumerate(variants[ntag]['samples']):
                            result = samplesi.copy()
                            g2 = v['samples'][i]
                            for j,samplej in enumerate(g2):
                                if g2[j] and g2[j]>0:
                                    result[j] = g2[j]
                            # print(i,samplesi,v['samples'][i],result)
                else:
                    variants[ntag] = v


        print("into")
        for key,v in variants.items():
            print(f"{key} AC={v['AC']}")

        self.assertEqual(len(variants),5)
        self.assertEqual(variants['10134532:G/T']['AC'],18)
        # Adjust TYPE field to set snp/mnp/ins/del
        key = None
        v = None
        for key,v in variants.items():
            ref_len = len(v['ref1'])
            aln_len = len(v['algn'])
            type = None
            size = None
            if aln_len < ref_len:
                type = 'del'
                size = ref_len - aln_len
            elif aln_len > ref_len:
                type = 'ins'
                size = aln_len - ref_len
            elif aln_len == ref_len:
                if ref_len == 1:
                    type = 'snp'
                else:
                    type = 'mnp'
                size = aln_len
            assert(size > 0)
            variants[key]['type'] = type
            variants[key]['size'] = size
            # Set origin
            print(v)
            variants[key]['origin'] = f'{var.name}:{var.pos}'
        # print(json.dumps(variants,indent=4))

        # handle deletions. If ref length is larger than the WF matched
        # allele length make this a missing genotype for all individual
        # SNP/MNP calls that match the allele index and fall inside the
        # deletion (region).
        #
        # The idea is that when a deletion exists for a sample there is
        # no way a SNP/MNP gets called in that sample. So for this deletion:
        # {'pos0': 10134514, 'ref0': 'GGAGAATCCCAATTGATGG', 'alt0': 'GG', 'ref1': 'GAGAATCCCAATTGATGG', 'algn': 'G', 'pos1': 10134515, 'altidx': 3, 'relpos': 1, 'AC': 3, 'AF': 0.0340909, 'AN': 88, 'is_rev': False, 'samples': [[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, None], [0, 0], [0, 0], [0, 0], [1, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [1, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [1, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0]], 'type': 'del'}
        # run through all other alleles and see if the genotype calls overlap. If
        # the SNP/MNP overlaps the deletion make it a NULL call. It means running
        # through all variants for every deletion. After discussion with Erik we
        # turn all haplotypes to NULL.
        type = None
        for key,v in variants.items():
            if v['type'] == 'del':
                # for every deletion
                del_ref_len = len(v['ref1'])
                del_aln_len = len(v['algn'])
                # del_len = del_ref_len - del_aln_len
                del_pos1 = v['pos1']
                del_size = v['size']
                del_start_pos = del_pos1 + del_aln_len
                # Make a range from the start of the deletion to the end
                check_range = range(del_start_pos, del_start_pos + del_size)
                check_samples = v['samples']
                for key2,v2 in variants.items():
                    if v2['type'] == 'snp' or v2['type'] == 'mnp':
                        # for alignment check all SNPs/MNPs
                        pos1 = v2['pos1']
                        pos2 = pos1 + v2['size']
                        if pos1 in check_range or pos2 in check_range:
                            # compare all genotypes
                            for i,sample in enumerate(v2['samples']):
                                del_sample = check_samples[i]
                                nullify = False
                                if 1 in del_sample and 1 in sample:
                                    nullify = True

                                    print(i,sample,del_sample,nullify)
                                    if nullify:
                                        # v2['samples'][i] = [None if item == 1 else item for item in sample]
                                        v2['samples'][i] = [None for item in sample]

        # Recompute AC and AF using the actual genotypes
        print("WIP")
        print(variants)

if __name__ == '__main__':
    unittest.main()
