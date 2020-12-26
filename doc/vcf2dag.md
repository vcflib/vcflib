% VCF2DAG(1) vcf2dag (vcflib) | vcf2dag (VCF unknown)
% Erik Garrison and vcflib contributors

# NAME

vcf2dag

# SYNOPSIS

usage: ./build/vcf2dag [options] [<vcf file>]

# DESCRIPTION

options: -r, --reference FILE FASTA reference file.

# OPTIONS

```


Modify the VCF file so that homozygous regions are included as REF/. calls.
For each ref and alt allele, assign an index.  These steps are sufficient to
enable use of the VCF as a DAG (specifically a partially-ordered graph).

```



# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# OTHER

## Source code

[vcf2dag.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcf2dag.cpp)

# LICENSE

Copyright 2011-2020 (C) Erik Garrison and vcflib contributors. MIT licensed.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
