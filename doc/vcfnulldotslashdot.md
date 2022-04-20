% VCFNULLDOTSLASHDOT(1) vcfnulldotslashdot | Convert VCF . to ./.
% Erik Garrison and other vcflib contributors

# NAME

vcfnulldotslashdot converts single dots to ./.

This is useful for some downstream analysis tools.

# SYNOPSIS

**vcf2tsv** < *file*

# DESCRIPTION

vcfnulldotslashdot converts single dots to ./.

This is useful for some downstream analysis tools.

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# EXAMPLES


<!--

    >>> from pytest.rtest import run_stdout, head, cat

-->

vcf2tsv converts a VCF to VCF

```python

>>> cat("../scripts/vcfnulldotslashdot data/issue_307_vcfnulldotslashdot.vcf")
##fileformat=VCFv4.1
##fileDate=20210209
##source=freeBayes v0.9.21
##reference=ahy.fa
##phasing=none
##filter="TYPE = snp & QUAL > 30 & AF > 0.05 & AF < 0.95 genotypes filtered with: GQ > 20"
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">
##INFO=<ID=LEN,Number=A,Type=Integer,Description="allele length">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	rab-field.AH13
contig71493_103113	50	.	G	A	79.0603	.	AC=5;AF=0.078125;LEN=1;TYPE=snp	GT	./.:0,0,0

```

## Source code

[vcfnulldotslashdot](../scripts/vcfnulldotslashdot)

# LICENSE

Copyright 2021 (C) Erik Garrison and vcflib contributors. MIT licensed.
