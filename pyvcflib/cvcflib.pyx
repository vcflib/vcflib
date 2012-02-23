# expose vcflib's VariantCallFile to Python as "VariantFile"
cdef class VariantFile:
    cdef VariantCallFile *vcffile_ptr
    
    def __init__(self, vcf_file):
        self.vcffile_ptr = new VariantCallFile()
        self.vcffile_ptr.openFile(string(vcf_file))
        self.fn = vcf_file