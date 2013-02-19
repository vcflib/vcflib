#!/usr/bin/env python

import sys
import time
from vcflib import PyVariantCallFile


if __name__ == '__main__':
    fn = sys.argv[1]
    vcf = PyVariantCallFile(fn)
    before = time.time()
    print len(vcf)
    after = time.time()
    print after-before
    print vcf.infoIds
    print vcf.formatIds
    print vcf.filterIds
    print vcf.infoTypes
    print vcf.formatTypes
    print vcf.infoCounts
    print vcf.formatCounts
    print vcf.parseSamples
    vcf = PyVariantCallFile(fn)
    vcf.parseSamples = False
    before = time.time()
    print len(vcf)
    after = time.time()
    print after-before
    vcf = PyVariantCallFile(fn)
    for variant in vcf:
        print variant