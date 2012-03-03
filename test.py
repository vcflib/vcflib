import pyvcflib
import sys

vcf_file = pyvcflib.VariantFile(sys.argv[1])
for var in vcf_file:
    # print var.chrom, var.pos, var.id, var.ref, \
    #       var.alt, var.qual, var.filter, var.info, \
    #       var.format, var.aaf, var.pi, var.num_alleles, \
    #       var.num_hom_ref, var.num_het, var.num_hom_alt, \
    #       var.num_unknown, var.num_valid, var.num_alt_alleles
          #var.samples \
    print var.chrom, var.pos, var.id, var.ref, \
        var.alt, var.num_hom_ref, \
        var.num_het, var.num_hom_alt, var.gts[0:10], var.gt_types[0:10]
        #var.num_unknown, var.num_valid, var.num_alt_alleles