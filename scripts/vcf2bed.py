#!/usr/bin/env python3
# Transform VCF to BED file
#
import sys

for line in sys.stdin:
    if line.startswith('#'):
        continue
    fields = line.strip().split()
    # VCF is 1-based, BED is 0-based half open
    # print out chrom, start, end, 
    chrom = fields[0]
    start = int(fields[1]) - 1
    span = len(fields[3]) # handle multi-base alleles
    end = start + span
    name = fields[2]
    print("\t".join(str(x) for x in [chrom, start, end, name]))
