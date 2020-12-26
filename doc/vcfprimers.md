% VCFPRIMERS(1) vcfprimers (vcflib) | vcfprimers (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcfprimers

# SYNOPSIS

usage: ./build/vcfprimers [options] <vcf file>

# DESCRIPTION

options: -f, --fasta-reference FASTA reference file to use to obtain primer sequences -l, --primer-length The length of the primer sequences on each side of the variant

# OPTIONS

```


For each VCF record, extract the flanking sequences, and write them to stdout as FASTA
records suitable for alignment.  This tool is intended for use in designing validation
experiments.  Primers extracted which would flank all of the alleles at multi-allelic
sites.  The name of the FASTA "reads" indicates the VCF record which they apply to.
The form is >CHROM_POS_LEFT for the 3' primer and >CHROM_POS_RIGHT for the 5' primer,
for example:

>20_233255_LEFT
CCATTGTATATATAGACCATAATTTCTTTATCCAATCATCTGTTGATGGA
>20_233255_RIGHT
ACTCAGTTGATTCCATACCTTTGCCATCATGAATCATGTTGTAATAAACA

```

# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcfprimers.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcfprimers.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
