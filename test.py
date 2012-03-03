import pyvcflib
import sys

vcf_file = pyvcflib.VariantFile(sys.argv[1])
for var in vcf_file:
    print var.chrom, var.pos, var.id, var.ref, \
          var.alt, var.qual, var.filter, var.info, \
          var.format, var.aaf, var.pi, var.samples
