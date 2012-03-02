import pyvcflib

vcf_file = pyvcflib.VariantFile("sample.vcf")
for var in vcf_file:
    print var.chrom, var.pos, var.id, var.ref, \
          var.alt, var.qual, var.filter, var.info, \
          var.format
    # 
    for sample in var.samples:
        print sample + ":" + str(var.samples[sample]['GT'][0])