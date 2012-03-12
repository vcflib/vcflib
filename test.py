import pyvcflib as pvcf
import sys

vcf_file = pvcf.VariantFile(sys.argv[1])

# echo the file
print vcf_file.header
for var in vcf_file:
    print var

# print individual attributes
for var in pvcf.VariantFile(sys.argv[1]):
    print "\t".join(str(s) for s in [var.chrom, var.pos, var.id, var.ref, \
        var.alt, var.qual, var.filter, var.info,\
        var.format, var.var_type, var.var_subtype,\
        var.num_hom_ref, var.num_het, var.num_hom_alt, \
        var.num_unknown, var.num_valid, var.num_alt_alleles, \
        var.gts[0:10], var.gt_types[0:10], var.gt_phases[0:10], \
        var.samples])
