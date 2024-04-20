% VCF2DAG(1) vcf2dag (vcflib) | vcf2dag (VCF transformation)
% Erik Garrison and vcflib contributors

# NAME

**vcf2dag**

# SYNOPSIS

**vcf2dag** [options] [<vcf file>]

# DESCRIPTION

Modify VCF to be able to build a directed acyclic graph (DAG)



# OPTIONS

```

options:
    -r, --reference FILE         FASTA reference file.

Modify the VCF file so that homozygous regions are included as REF/. calls.
For each ref and alt allele, assign an index.  These steps are sufficient to
enable use of the VCF as a DAG (specifically a partially-ordered graph).

Type: transformation

```





# EXIT VALUES

**0**
: Success

**not 0**
: Failure

# SEE ALSO



[vcflib](./vcflib.md)(1)



# OTHER

## Source code

[vcf2dag.cpp](https://github.com/vcflib/vcflib/blob/master/src/vcf2dag.cpp)

# LICENSE

Copyright 2011-2024 (C) Erik Garrison and vcflib contributors. MIT licensed.
Copyright 2020-2024 (C) Pjotr Prins.

<!--
  Created with ./scripts/bin2md.rb scripts/bin2md-template.erb
-->
